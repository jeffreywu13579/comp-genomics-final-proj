#!/usr/bin/env python3
import os
import sys
import time
import argparse
import requests
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def query_cbioportal(path, params=None, method='GET', json_body=None):
    base = "https://www.cbioportal.org/api"
    url = base + path
    headers = {"Accept": "application/json"}
        
    if method == 'GET':
        r = requests.get(url, headers=headers, params=params)
    else:
        r = requests.post(url, headers=headers, json=json_body)
    
    try:
        r.raise_for_status()
        response_data = r.json()
        return response_data
    except requests.exceptions.RequestException as e:
        raise


def fetch_mutation_data(mut_profile, sample_list_id, entrez_ids, symbols):
    mut_payload = {
        "sampleListId": sample_list_id,
        "entrezGeneIds": entrez_ids
    }   
    try:
        mut_data = query_cbioportal(
            f"/molecular-profiles/{mut_profile}/mutations/fetch", 
            method='POST', 
            json_body=mut_payload
        )
        return mut_data
    except Exception as e:
        return fetch_mutations_individually(mut_profile, sample_list_id, entrez_ids, symbols)

def fetch_mutations_individually(mut_profile, sample_list_id, entrez_ids, symbols):
    mut_records = []
    for gid in entrez_ids:
        try:
            data = query_cbioportal(
                f"/molecular-profiles/{mut_profile}/mutations",
                params={"entrezGeneId": gid, "sampleListId": sample_list_id}
            )
            mut_records.extend(data)
            time.sleep(0.5)
        except Exception as e:
            pass
    return mut_records


def fetch_cbio_data(study_id, max_genes=50):
    prefix = "_".join(study_id.split("_")[:2])
    sample_lists = query_cbioportal("/sample-lists", params={"studyId": study_id})
    all_sample_list_ids = [sl['sampleListId'] for sl in sample_lists]
    sample_list_id = next((sl for sl in all_sample_list_ids if sl.startswith(prefix) and sl.endswith("_all")), None)
    if not sample_list_id:
        raise ValueError("No valid '_all' sample list found for this study")

    try:
        mut_sig_genes = query_cbioportal(f"/studies/{study_id}/significantly-mutated-genes")
        mut_sig_genes.sort(key=lambda x: x.get('qValue', 1.0))
        top_genes = [g['hugoGeneSymbol'] for g in mut_sig_genes[:max_genes]]
        if not top_genes:
            top_genes = [
                "EGFL7", "GATA6", "ASXL2", "RALGDS", "ZNF24", "H3C13", "PDGFB", "HOOK3",
                "SFRP2", "PTPRC", "GEN1", "ZRSR2", "STAG2", "IGH", "MRE11", "H3C2",
                "HOXA13", "STAT4", "CREB3L1", "SOX17", "GPS2", "MAFB", "PALB2", "MN1",
                "MLF1", "DIS3", "WDR90", "FOLR1", "SESN3", "DNM2", "CLTCL1", "DIS3L2",
                "AFF1", "PTPRO", "ELF3", "TCEA1", "ATF1", "CRLF2", "CTDNEP1", "MLH1"
            ]
    except Exception as e:
        top_genes = [
                "EGFL7", "GATA6", "ASXL2", "RALGDS", "ZNF24", "H3C13", "PDGFB", "HOOK3",
                "SFRP2", "PTPRC", "GEN1", "ZRSR2", "STAG2", "IGH", "MRE11", "H3C2",
                "HOXA13", "STAT4", "CREB3L1", "SOX17", "GPS2", "MAFB", "PALB2", "MN1",
                "MLF1", "DIS3", "WDR90", "FOLR1", "SESN3", "DNM2", "CLTCL1", "DIS3L2",
                "AFF1", "PTPRO", "ELF3", "TCEA1", "ATF1", "CRLF2", "CTDNEP1", "MLH1"
            ]
    top_genes = top_genes[:max_genes]

    genes = []
    page = 0
    while True:
        res = query_cbioportal("/genes", params={"pageSize": 1000, "pageNumber": page})
        if not res:
            break
        genes.extend(res)
        page += 1
    genes_df = pd.DataFrame(genes)
    
    genes_df = genes_df[genes_df['hugoGeneSymbol'].isin(top_genes)]
    entrez_ids = genes_df['entrezGeneId'].dropna().astype(int).tolist()
    symbols = genes_df.set_index('entrezGeneId')['hugoGeneSymbol'].to_dict()
    
    expr_profile = f"{prefix}_rna_seq_v2_mrna"
    expr_records = []
    for i, gid in enumerate(entrez_ids):
        try:
            data = query_cbioportal(
                f"/molecular-profiles/{expr_profile}/molecular-data",
                params={"entrezGeneId": gid, "sampleListId": sample_list_id}
            )
            if not data:
                continue
            expr_records.extend(data)
            time.sleep(0.2)
        except Exception as e:
            pass
    expr_df = pd.DataFrame(expr_records)
    assert not expr_df.empty, "No expression data was retrieved"
    
    required_cols = {'sampleId', 'entrezGeneId', 'value'}
    missing_cols = required_cols - set(expr_df.columns)
    if missing_cols:
        raise ValueError(f"Expression data missing columns: {missing_cols}")
    
    expr_df = expr_df.dropna(subset=['sampleId', 'entrezGeneId', 'value'])
    expr_mat = expr_df.pivot(index='sampleId', columns='entrezGeneId', values='value')
    expr_mat = expr_mat.reindex(columns=entrez_ids).fillna(0)
    print(f"Expression matrix shape: {expr_mat.shape}")

    cna_profile = f"{prefix}_gistic"
    cna_data = []
    for i, gid in enumerate(entrez_ids):
        try:
            data = query_cbioportal(
                f"/molecular-profiles/{cna_profile}/molecular-data",
                params={"entrezGeneId": gid, "sampleListId": sample_list_id}
            )
            if not data:
                continue
            cna_data.extend(data)
            time.sleep(0.2)
        except Exception as e:
            pass
    
    cna_df = pd.DataFrame(cna_data)
    cna_mat = cna_df.pivot(index='sampleId', columns='entrezGeneId', values='value').reindex(columns=entrez_ids).fillna(0)
    amp_mat = (cna_mat == 2).astype(int)
    del_mat = (cna_mat == -2).astype(int)

    mut_profile = f"{prefix}_mutations"
    mut_records = fetch_mutation_data(mut_profile, sample_list_id, entrez_ids, symbols)
    
    samples = sorted(set(expr_mat.index) & set(cna_mat.index))
    mis_mat = pd.DataFrame(0, index=samples, columns=entrez_ids)
    non_mat = pd.DataFrame(0, index=samples, columns=entrez_ids)
    
    mutation_counts = {gid: 0 for gid in entrez_ids}
    if mut_records:
        for mut in mut_records:
            sid = mut.get('sampleId')
            gid = mut.get('entrezGeneId')
            
            if sid is None or gid is None:
                continue
                
            if sid in mis_mat.index and gid in mis_mat.columns:
                mtype = mut.get('mutationType', '') 
                if not mtype:
                    mtype = mut.get('mutationClassification', '')
                if 'Missense' in mtype:
                    mis_mat.at[sid, gid] = 1
                    mutation_counts[gid] = mutation_counts.get(gid, 0) + 1
                elif 'Nonsense' in mtype:
                    non_mat.at[sid, gid] = 1
                    mutation_counts[gid] = mutation_counts.get(gid, 0) + 1

    all_samples = sorted(set(expr_mat.index) & set(amp_mat.index) & set(del_mat.index) & set(mis_mat.index) & set(non_mat.index))
    
    expr_mat = expr_mat.loc[all_samples]
    amp_mat = amp_mat.loc[all_samples]
    del_mat = del_mat.loc[all_samples]
    mis_mat = mis_mat.loc[all_samples]
    non_mat = non_mat.loc[all_samples]
    
    return expr_mat, amp_mat, del_mat, mis_mat, non_mat, entrez_ids, symbols, all_samples


def run_notears(expr, amp, dele, mis, non, entrez_ids, symbols, num_edges_drawn=30):
    from notears.linear import notears_linear
    X = np.hstack([expr.values, amp.values, dele.values, mis.values, non.values])
    feature_names = []
    for gid in entrez_ids:
        sym = symbols.get(gid, str(gid))
        feature_names.extend([f"{sym}_expr", f"{sym}_amp", f"{sym}_del", f"{sym}_mis", f"{sym}_non"])

    print("Running NOTEARS")
    W = notears_linear(X, lambda1=0.01, loss_type='l2')
    os.makedirs("results", exist_ok=True)
    pd.DataFrame(W, index=feature_names, columns=feature_names).to_csv("results/adjacency_matrix.csv")

    max_abs_weight = np.max(np.abs(W[W!=0]))
    if max_abs_weight > 0:
        W_norm = W / max_abs_weight 
    else:
        W_norm = W

    edges = [(feature_names[i], feature_names[j], W_norm[i, j], abs(W_norm[i, j]))
             for i in range(len(feature_names)) for j in range(len(feature_names)) if i != j and W_norm[i, j] != 0]
    top_edges = sorted(edges, key=lambda x: x[3], reverse=True)[:50]
    pd.DataFrame(top_edges, columns=['source', 'target', 'weight', 'abs_weight']).to_csv("results/top_edges.csv", index=False)

    G = nx.DiGraph()
    for src, dst, w, _ in top_edges[:num_edges_drawn]:
        G.add_edge(src, dst, weight=w)
    
    pos = nx.spring_layout(G, k=0.5, seed=42)
    plt.figure(figsize=(14, 14)) 
    
    expr_nodes = [n for n in G.nodes() if n.endswith('_expr')]
    amp_nodes = [n for n in G.nodes() if n.endswith('_amp')]
    del_nodes = [n for n in G.nodes() if n.endswith('_del')]
    mis_nodes = [n for n in G.nodes() if n.endswith('_mis')]
    non_nodes = [n for n in G.nodes() if n.endswith('_non')]
    
    nx.draw_networkx_nodes(G, pos, nodelist=expr_nodes, node_size=500, node_color='lightgreen', alpha=0.8)
    nx.draw_networkx_nodes(G, pos, nodelist=amp_nodes, node_size=500, node_color='lightcoral', alpha=0.8)
    nx.draw_networkx_nodes(G, pos, nodelist=del_nodes, node_size=500, node_color='lightblue', alpha=0.8)
    nx.draw_networkx_nodes(G, pos, nodelist=mis_nodes, node_size=500, node_color='gold', alpha=0.8)
    nx.draw_networkx_nodes(G, pos, nodelist=non_nodes, node_size=500, node_color='mediumpurple', alpha=0.8)
    
    for node_list in [expr_nodes, amp_nodes, del_nodes, mis_nodes, non_nodes]:
        nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_size=500, 
                              node_color='none', edgecolors='black', linewidths=1)
    
    nx.draw_networkx_labels(G, pos, font_size=8)
    
    positive_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] > 0]
    negative_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] < 0]
    
    pos_widths = [2.5 * abs(G[u][v]['weight']) for u, v in positive_edges]
    neg_widths = [2.5 * abs(G[u][v]['weight']) for u, v in negative_edges]
    
    if positive_edges:
        nx.draw_networkx_edges(
            G, pos, 
            edgelist=positive_edges, 
            edge_color='green', 
            arrows=True,
            width=pos_widths,
            arrowsize=15
        )
    
    if negative_edges:  
        nx.draw_networkx_edges(
            G, pos, 
            edgelist=negative_edges, 
            edge_color='red', 
            arrows=True,
            width=neg_widths,
            arrowsize=15
        )
    
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='green', lw=2, label='Positive Influence'),
        Line2D([0], [0], color='red', lw=2, label='Negative Influence'),
        Line2D([0], [0], marker='o', color='lightgreen', label='Expression', markerfacecolor='lightgreen', markersize=10),
        Line2D([0], [0], marker='o', color='lightcoral', label='Amplification', markerfacecolor='lightcoral', markersize=10),
        Line2D([0], [0], marker='o', color='lightblue', label='Deletion', markerfacecolor='lightblue', markersize=10),
        Line2D([0], [0], marker='o', color='gold', label='Missense Mutation', markerfacecolor='gold', markersize=10),
        Line2D([0], [0], marker='o', color='mediumpurple', label='Nonsense Mutation', markerfacecolor='mediumpurple', markersize=10)
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.title(f"Gene Interactions Network (Top {len(G.edges())} edges)", fontsize=16)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig("results/causal_network.png", dpi=300)
    print("Saved visualization to 'results/causal_network.png'")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("study", type=str, help="Study name: BRCA, COADREAD, BLCA")
    parser.add_argument("--max-genes", type=int, default=25, help="Max number of top genes to include (default: 25)")
    parser.add_argument("--num-edges-drawn", type=int, default=30, help="Max number of edges drawn (default: 30)")
    args = parser.parse_args()
    
    study_map = {
        "BRCA": "brca_tcga_pan_can_atlas_2018",
        "COADREAD": "coadread_tcga_pan_can_atlas_2018",
        "BLCA": "blca_tcga_pan_can_atlas_2018"
    }
    study_id = study_map.get(args.study.upper())
    
    print(f"Analyzing {study_id} with max {args.max_genes} genes")
    try:
        expr, amp, dele, mis, non, entrez_ids, symbols, samples = fetch_cbio_data(study_id, args.max_genes)
        assert len(expr) != 0, "No valid expression data found."
        print("Data retrieval complete.")
        run_notears(expr, amp, dele, mis, non, entrez_ids, symbols, args.num_edges_drawn)
    except Exception as e:
        print(f"Error encountered: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()