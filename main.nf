#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * functionalhomology
 * E. coli K-12 vs B. subtilis 168: ESM-2 embedding similarity versus DIAMOND
 * sequence similarity, evaluated with experimental GO annotations.
 */

process DIAMOND {
    tag "${proteome_a.baseName} vs ${proteome_b.baseName}"
    label 'process_medium'
    publishDir "${params.outdir}/diamond", mode: 'copy'

    input:
    path proteome_a
    path proteome_b

    output:
    path 'diamond_hits.tsv', emit: hits

    script:
    """
    set -euo pipefail
    diamond makedb --in "${proteome_b}" -d target_db --threads ${task.cpus}
    diamond blastp \\
        -q "${proteome_a}" \\
        -d target_db \\
        --very-sensitive \\
        -e ${params.evalue_cutoff} \\
        --max-target-seqs 0 \\
        --max-hsps 1 \\
        --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident evalue bitscore \\
        -o diamond_hits.tsv
    echo "DIAMOND reported \$(wc -l < diamond_hits.tsv) HSPs"
    """
}

process EMBED {
    tag "${name}"
    label 'process_embeddings'
    publishDir "${params.outdir}/embeddings", mode: 'copy'

    input:
    tuple path(fasta), val(name)

    output:
    tuple val(name), path("${name}.npy"), path("${name}.ids.txt"), path("${name}.meta.json"), emit: embeddings

    script:
    """
    set -euo pipefail
    embed.py \\
        --fasta "${fasta}" \\
        --outdir . \\
        --name "${name}" \\
        --model "${params.esm_model}" \\
        --batch-size ${params.embed_batch_size} \\
        --max-residues ${params.max_residues}
    """
}

process FIND_PAIRS {
    tag "top-k=${params.top_k} cos>=${params.cosine_threshold}"
    label 'process_low'
    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    tuple path(emb_a), path(ids_a)
    tuple path(emb_b), path(ids_b)
    path diamond_hits

    output:
    path 'pairs.tsv', emit: pairs
    path 'diamond_pairs.tsv', emit: diamond
    path 'random_pairs.tsv', emit: random
    path 'pairs_summary.json', emit: summary

    script:
    """
    set -euo pipefail
    find_pairs.py \\
        --embeddings-a "${emb_a}" \\
        --ids-a "${ids_a}" \\
        --embeddings-b "${emb_b}" \\
        --ids-b "${ids_b}" \\
        --diamond-hits "${diamond_hits}" \\
        --outdir . \\
        --top-k ${params.top_k} \\
        --cosine-threshold ${params.cosine_threshold} \\
        --evalue-cutoff ${params.evalue_cutoff} \\
        --seed ${params.seed}
    """
}

process GO_EVAL {
    tag "experimental GO"
    label 'process_low'
    publishDir "${params.outdir}/go_eval", mode: 'copy'

    input:
    path pairs
    path diamond_pairs
    path random_pairs
    path go_ann
    path go_obo

    output:
    path 'pairs_scored.tsv', emit: scored
    path 'eval_summary.json', emit: summary

    script:
    """
    set -euo pipefail
    go_eval.py \\
        --pairs "${pairs}" \\
        --diamond "${diamond_pairs}" \\
        --random "${random_pairs}" \\
        --go-ann "${go_ann}" \\
        --go-obo "${go_obo}" \\
        --outdir . \\
        --cosine-threshold ${params.cosine_threshold}
    """
}

process REPORT {
    tag "HTML report"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path summary
    path scored
    path pairs
    path fasta_a
    path fasta_b

    output:
    path 'report.html'

    script:
    """
    set -euo pipefail
    make_report.py \\
        --summary "${summary}" \\
        --pairs-scored "${scored}" \\
        --pairs "${pairs}" \\
        --fasta-a "${fasta_a}" \\
        --fasta-b "${fasta_b}" \\
        --out report.html \\
        --model "${params.esm_model}" \\
        --top-k ${params.top_k} \\
        --cosine-threshold ${params.cosine_threshold} \\
        --seed ${params.seed}
    """
}

workflow {
log.info """
============================================================
 functionalhomology
------------------------------------------------------------
 proteome A     : ${params.proteome_a}
 proteome B     : ${params.proteome_b}
 GO annotations : ${params.go_ann}
 GO ontology    : ${params.go_obo}
 ESM model      : ${params.esm_model}
 top-k          : ${params.top_k}
 cosine cutoff  : ${params.cosine_threshold}
 seed           : ${params.seed}
 output         : ${params.outdir}
============================================================
"""

if (!params.proteome_a || !params.proteome_b || !params.go_ann || !params.go_obo) {
    error "Missing required input. Set --proteome_a, --proteome_b, --go_ann and --go_obo " +
          "(or use -profile test)."
}

    proteome_a = file(params.proteome_a, checkIfExists: true)
    proteome_b = file(params.proteome_b, checkIfExists: true)
    go_ann     = file(params.go_ann,     checkIfExists: true)
    go_obo     = file(params.go_obo,     checkIfExists: true)

    DIAMOND(proteome_a, proteome_b)

    // One EMBED task per proteome, executed in parallel.
    EMBED(Channel.of(
        tuple(proteome_a, 'proteome_a'),
        tuple(proteome_b, 'proteome_b')
    ))
    embeddings_a = EMBED.out.embeddings.filter { name, npy, ids, meta -> name == 'proteome_a' }
    embeddings_b = EMBED.out.embeddings.filter { name, npy, ids, meta -> name == 'proteome_b' }

    FIND_PAIRS(
        embeddings_a.map { name, npy, ids, meta -> tuple(npy, ids) },
        embeddings_b.map { name, npy, ids, meta -> tuple(npy, ids) },
        DIAMOND.out.hits
    )

    GO_EVAL(
        FIND_PAIRS.out.pairs,
        FIND_PAIRS.out.diamond,
        FIND_PAIRS.out.random,
        go_ann,
        go_obo
    )

    REPORT(
        GO_EVAL.out.summary,
        GO_EVAL.out.scored,
        FIND_PAIRS.out.pairs,
        proteome_a,
        proteome_b
    )
}
