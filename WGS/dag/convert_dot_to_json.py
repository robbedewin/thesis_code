import pydot
import json

dot_data = """
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans, fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
    0[label = "all", color = "0.24 0.6 0.85", style="rounded"];
    1[label = "multiqc", color = "0.51 0.6 0.85", style="rounded"];
    2[label = "fastqc_raw", color = "0.11 0.6 0.85", style="rounded"];
    3[label = "fastqc_trimmed", color = "0.04 0.6 0.85", style="rounded"];
    4[label = "trim_reads_for_fastqc", color = "0.61 0.6 0.85", style="rounded"];
    5[label = "samtools_idxstats", color = "0.47 0.6 0.85", style="rounded"];
    6[label = "apply_bqsr", color = "0.15 0.6 0.85", style="rounded"];
    7[label = "mark_duplicates", color = "0.37 0.6 0.85", style="rounded"];
    8[label = "samtools_sort", color = "0.12 0.6 0.85", style="rounded"];
    9[label = "bwa_map", color = "0.64 0.6 0.85", style="rounded"];
    10[label = "trim_reads", color = "0.63 0.6 0.85", style="rounded"];
    11[label = "get_reference", color = "0.52 0.6 0.85", style="rounded"];
    12[label = "bwa_index", color = "0.20 0.6 0.85", style="rounded"];
    13[label = "base_recal_table", color = "0.39 0.6 0.85", style="rounded"];
    14[label = "genome_faidx", color = "0.49 0.6 0.85", style="rounded"];
    15[label = "genome_dict", color = "0.00 0.6 0.85", style="rounded"];
    16[label = "get_common_dbSNP", color = "0.65 0.6 0.85", style="rounded"];
    17[label = "samtools_index", color = "0.03 0.6 0.85", style="rounded"];
    18[label = "samtools_stats", color = "0.41 0.6 0.85", style="rounded"];
    19[label = "ascat", color = "0.59 0.6 0.85", style="rounded"];
    20[label = "mutect2", color = "0.23 0.6 0.85", style="rounded"];
    21[label = "mutect_pon3", color = "0.07 0.6 0.85", style="rounded"];
    22[label = "mutect_pon2", color = "0.40 0.6 0.85", style="rounded"];
    23[label = "create_sample_map", color = "0.31 0.6 0.85", style="rounded"];
    24[label = "mutect_pon1", color = "0.28 0.6 0.85", style="rounded"];
    25[label = "select_pass_variants", color = "0.13 0.6 0.85", style="rounded"];
    26[label = "filter_mutect2", color = "0.17 0.6 0.85", style="rounded"];
    27[label = "calculate_contamination", color = "0.32 0.6 0.85", style="rounded"];
    28[label = "get_pileup_summaries", color = "0.16 0.6 0.85", style="rounded"];
    29[label = "germline_resource_chr", color = "0.25 0.6 0.85", style="rounded"];
    30[label = "germline_resource", color = "0.53 0.6 0.85", style="rounded"];
    31[label = "learn_read_orientation_model", color = "0.08 0.6 0.85", style="rounded"];
    32[label = "liftover_chm13_to_grch38_vcf", color = "0.19 0.6 0.85", style="rounded"];
    33[label = "funcotate_variants", color = "0.36 0.6 0.85", style="rounded"];
    34[label = "call_somatic_structural_variants", color = "0.48 0.6 0.85", style="rounded"];
    35[label = "gridss_somatic_filter", color = "0.09 0.6 0.85", style="rounded"];
    36[label = "generate_pon", color = "0.60 0.6 0.85", style="rounded"];
    37[label = "gridss_normal", color = "0.55 0.6 0.85", style="rounded"];
    38[label = "get_exclude_regions", color = "0.43 0.6 0.85", style="rounded"];
    39[label = "svaba", color = "0.29 0.6 0.85", style="rounded"];
    17 -> 0
    32 -> 0
    34 -> 0
    33 -> 0
    35 -> 0
    20 -> 0
    1 -> 0
    25 -> 0
    19 -> 0
    39 -> 0
    5 -> 1
    3 -> 1
    2 -> 1
    7 -> 1
    18 -> 1
    4 -> 3
    17 -> 5
    6 -> 5
    7 -> 6
    11 -> 6
    13 -> 6
    8 -> 7
    9 -> 8
    11 -> 8
    10 -> 9
    12 -> 9
    11 -> 9
    11 -> 12
    15 -> 13
    14 -> 13
    7 -> 13
    11 -> 13
    16 -> 13
    11 -> 14
    11 -> 15
    6 -> 17
    6 -> 18
    6 -> 19
    6 -> 20
    21 -> 20
    11 -> 21
    22 -> 21
    23 -> 22
    11 -> 22
    24 -> 23
    11 -> 24
    6 -> 24
    26 -> 25
    27 -> 26
    31 -> 26
    20 -> 26
    28 -> 27
    6 -> 28
    29 -> 28
    30 -> 29
    20 -> 31
    25 -> 32
    32 -> 33
    6 -> 34
    34 -> 35
    36 -> 35
    37 -> 36
    6 -> 37
    38 -> 37
    6 -> 39
}
"""

# Parse the dot data
(graph,) = pydot.graph_from_dot_data(dot_data)

nodes = []
links = []

for node in graph.get_nodes():
    nodes.append({
        "id": node.get_name(),
        "label": node.get_label().strip('"') if node.get_label() else node.get_name()
    })

for edge in graph.get_edges():
    links.append({
        "source": edge.get_source(),
        "target": edge.get_destination()
    })

# Create the JSON structure
dag_json = {
    "nodes": nodes,
    "links": links
}

# Save to a JSON file
with open('dag/dag_py_converted.json', 'w') as f:
    json.dump(dag_json, f, indent=4)
