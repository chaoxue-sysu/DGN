# -*- coding: utf-8 -*-
# @Author  : Xue Chao
# @Time    : 2024/01/11 11:37
# @Function: Main module of DGN (disease-associated gene network module).
#
import argparse
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd
import scipy
import seaborn
from CSCORE import CSCORE_IRLS
from gprofiler import GProfiler
from matplotlib import pyplot as plt, font_manager
from pyecharts.charts import Graph
from pyecharts import options as opts
from scipy import stats
from scipy.optimize import fsolve
from scipy.stats import zscore

from util import index_in_gene_symbol_of_gene_list, log, cpm, make_dir, cal_hist_density, FONTS_DIR, kggsee_dese, \
    run_command, KGGSEE, pvalue_adjust, need_trans_gene_name, kggsee_rez, cell_corr_heatmap, \
    cluster_df, heatmap_color_rstyle, heatmap_color_rstyle_single, jaccard, go_enrich_plot, GeneNetwork, KGGSEE_DIR

import warnings
warnings.filterwarnings("ignore")
font_files = font_manager.findSystemFonts(fontpaths=FONTS_DIR)
for file in font_files:
    font_manager.fontManager.addfont(file)
plt.rcParams["font.sans-serif"] = "Arial"

PKG_NAME='DGN'
VERSION='1.0.0'

class DGN:
    '''
    DGN class
    '''
    def __init__(self,input_dict:{}):
        self.para=input_dict
        self.para['kggsee_dir']=KGGSEE_DIR
        out_dir=self.para['output_dir']
        self.node_score_dir = f'{out_dir}/node_score'# save: degree, mean expression.
        self.intermediate_dir = f'{out_dir}/intermediate'# save: k factor, resample corr coef.
        self.result_dir=f'{out_dir}/result'
        self.saved_ref_dir = f'{out_dir}/saved_ref'
        self.unified_genes_path=f'{self.intermediate_dir}/unified_genes.txt.gz'
        make_dir(self.node_score_dir,self.intermediate_dir)
        self.cell_k = {}
        self.module_p_cutoff = 0.1
        self.node_score_file_name='degree.k_adjusted.txt.gz'

    def __print_effect_para(self):
        msg='Effective options:'
        for k,v in self.para.items():
            if v is not None:
                msg+=f'\n --{k} {v}'
        print(msg)

    def extract_expr_top(expr: np.ndarray, gene_id: dict, n_top_genes: int, expr_norm_method: str, include_genes=None):
        all_genes = expr.shape[0]
        if all_genes < n_top_genes:
            n_top_genes = all_genes
        method_map = {'no': lambda x: x, 'cpm': cpm}
        norm_expr = method_map[expr_norm_method](expr)
        mean_expr = np.nanmean(norm_expr, axis=1)
        sorted_indices = np.argsort(mean_expr)
        top_indices = sorted_indices[-n_top_genes:][::-1]
        if include_genes is not None:
            include_idx = []
            genes = gene_id['gene_symbol']
            for i in range(len(genes)):
                if genes[i] in include_genes:
                    include_idx.append(i)
            remain_idx = np.intersect1d(top_indices, include_idx)
        else:
            remain_idx = top_indices
        fine_expr = expr[remain_idx, :]
        fine_gene_id = {}
        for k in gene_id.keys():
            fine_gene_id[k] = gene_id[k][remain_idx]
        return fine_expr, fine_gene_id

    def __select_sub_cell(self, cell_ids, max_cell):
        if len(cell_ids)<max_cell:
            max_cell=len(cell_ids)
        np.random.seed(10086)
        sub=np.random.choice(cell_ids,max_cell,replace=False)
        return sub

    def __read_expr_file(self,include_cells=None):
        exclude_cells=None
        na_fill_val = 0
        expr_file=self.para["expression_path"]
        cell_annot_file=self.para["cell_annot_path"]
        trans_gene_symbol = self.para['trans_gene_symbol']
        min_cell=self.para["min_cell"]
        max_cell = self.para["max_cell"]
        expr_df=pd.read_table(expr_file,header=0,index_col=0)
        expr_df = expr_df.loc[~expr_df.index.duplicated(), :]
        anno_df=pd.read_table(cell_annot_file,header=None,index_col=None)
        cell_types=sorted(anno_df[1].unique())
        if os.path.isfile(self.unified_genes_path):
            ug_df=pd.read_table(self.unified_genes_path)
            gene_id={c:ug_df[c].values for c in ug_df.columns}
        else:
            unique_genes=expr_df.index.tolist()
            if not trans_gene_symbol:
                if need_trans_gene_name(unique_genes[0]):
                    trans_gene_symbol=True
            if trans_gene_symbol:
                gene_idx, gene_syms = index_in_gene_symbol_of_gene_list(unique_genes)
                fine_genes = np.array(unique_genes)[gene_idx]
                fine_gene_sym = np.array(gene_syms)[gene_idx]
                log(f'raw {len(unique_genes)} genes; remain {len(fine_genes)} genes.')
            else:
                fine_genes = np.array(unique_genes)
                fine_gene_sym = np.copy(fine_genes)
            gene_id = {'gene_id': fine_genes, 'gene_symbol': fine_gene_sym}
            ug_df = pd.DataFrame(gene_id)
            ug_df.to_csv(self.unified_genes_path, sep='\t', index=False)
            unique_cells=set(expr_df.columns.tolist())
        gene_raw_id=gene_id['gene_id']
        log(f'loaded {len(gene_id["gene_symbol"])} unified genes.')
        for cell_name in cell_types:
            if include_cells is not None and cell_name not in include_cells:
                continue
            if exclude_cells is not None and cell_name in exclude_cells:
                continue
            cell_ids=anno_df.loc[anno_df[1]==cell_name,0]
            cell_ids=sorted(unique_cells.intersection(cell_ids))
            if len(cell_ids)<min_cell:
                continue
            if len(cell_ids)>max_cell:
                cell_ids=self.__select_sub_cell(cell_ids,max_cell)
            df=expr_df.loc[:,cell_ids]
            df = df.loc[gene_raw_id, :].fillna(na_fill_val)
            expr=df.values
            yield expr, re.sub('\s+','',cell_name), gene_id

    def __read_expr_dir(self,include_cells=None):
        exclude_cells=None
        expr_dir=self.para["expression_dir"]
        min_cell=self.para["min_cell"]
        max_cell = self.para["max_cell"]
        trans_gene_symbol = self.para['trans_gene_symbol']
        na_fill_val = 0
        expr_files=sorted(os.listdir(expr_dir))
        if os.path.isfile(self.unified_genes_path):
            ug_df=pd.read_table(self.unified_genes_path)
            gene_id={c:ug_df[c].values for c in ug_df.columns}
        else:
            unique_genes = set()
            for f in expr_files:
                try:
                    df = pd.read_table(f'{expr_dir}/{f}', header=0, index_col=0)
                except:
                    raise Exception(f'error unzip in {f}')
                if df.shape[1] < min_cell:
                    continue
                unique_genes = unique_genes.union(df.index.values)
            unique_genes = sorted(unique_genes)
            if not trans_gene_symbol:
                if need_trans_gene_name(unique_genes[0]):
                    trans_gene_symbol=True
            if trans_gene_symbol:
                gene_idx, gene_syms = index_in_gene_symbol_of_gene_list(unique_genes)
                fine_genes = np.array(unique_genes)[gene_idx]
                fine_gene_sym = np.array(gene_syms)[gene_idx]
                log(f'raw {len(unique_genes)} genes; remain {len(fine_genes)} genes.')
            else:
                fine_genes = np.array(unique_genes)
                fine_gene_sym = np.copy(fine_genes)
            gene_id = {'gene_id': fine_genes, 'gene_symbol': fine_gene_sym}
            ug_df=pd.DataFrame(gene_id)
            ug_df.to_csv(self.unified_genes_path,sep='\t',index=False)
        gene_raw_id=gene_id['gene_id']
        log(f'loaded {len(gene_id["gene_symbol"])} unified genes.')
        for f in expr_files:
            last_idx=1
            if f.endswith('.txt.gz'):
                last_idx=2
            arr=f.split('.')
            if len(arr)<2:
                last_idx=0
            cell_name = '.'.join(arr[:-last_idx])
            if include_cells is not None and cell_name not in include_cells:
                continue
            if exclude_cells is not None and cell_name in exclude_cells:
                continue
            df = pd.read_table(f'{expr_dir}/{f}', header=0, index_col=0)
            cell_ids=df.columns.values
            if len(cell_ids) < min_cell:
                continue
            if len(cell_ids) > max_cell:
                cell_ids=self.__select_sub_cell(cell_ids,max_cell)
            df = df.loc[~df.index.duplicated(), cell_ids]
            add_genes = set(gene_raw_id).union(df.index) - set(df.index)
            for g in add_genes:
                df.loc[g, :] = na_fill_val
            df = df.loc[gene_raw_id, :].fillna(na_fill_val)
            expr=df.values
            yield expr, re.sub('\s+','',cell_name), gene_id

    def __normalize_expr(self,expr):
        method_map={'no':lambda x:x,'cpm':cpm}
        normalize_expr=self.para["normalize_expr"]
        return method_map[normalize_expr](expr)

    def corr_matrix(self,expr,seq_depth):
        def cs_core(expr,depth):
            counts=expr.T
            net = CSCORE_IRLS(counts, depth)
            return net[0]
        method=self.para['edge_method']
        min_expr_val=self.para["min_expr_value"]
        zero_idx=np.where(np.nanmean(expr,axis=1)<min_expr_val)[0]
        ## The genes with extremely low expression levels lack credibility
        if min_expr_val>0:
            log(f'{len(zero_idx)} (of {expr.shape[0]}) genes with low expression values (mean<{min_expr_val}).')
        # if method!='cs-core':
        #     expr=self.__normalize_expr(expr)
        #     if self.para['log_trans_expr']:
        #         expr=np.log2(expr+1)
        if method=='cs-core':
            corr=cs_core(expr,seq_depth)
        if method=='pearson':
            corr=np.corrcoef(expr)
        if len(zero_idx)>0:
            corr[zero_idx,:]=0
            corr[:,zero_idx]=0
        np.fill_diagonal(corr, 1)
        return corr


    def __cal_k_factor(self,resample_corr,norm_range=[0, 1]):
        def k_equation(k, X, v):
            return np.nanmean(np.power(X, k)) - v

        datas = np.abs(resample_corr)
        st_idx, ed_idx = (np.array(norm_range) * datas.shape[0]).astype(int)
        sort_datas = np.sort(datas, axis=0)
        norm_datas = sort_datas[st_idx:ed_idx, :]
        v = np.nanmean(np.nanmean(norm_datas, axis=1))
        ks = []
        for i in range(norm_datas.shape[1]):
            k = fsolve(k_equation, 0.1, (norm_datas[:, i], v))
            ks.append(k[0])
        ks=np.array(ks)
        return ks

    def __plot_distribution(self,datas,cells,fig_path,info=None):
        sort_datas = np.sort(datas, axis=0)
        all_bg=sort_datas.flatten()
        nrow, ncol = datas.shape
        need_axes_num=ncol
        sub_axes_unit=3
        n = int(np.sqrt(need_axes_num) - 0.0001) + 1
        plot_ncol = n
        plot_nrow = int(need_axes_num/plot_ncol-0.0001)+1
        fig, axes = plt.subplots(plot_nrow, plot_ncol, figsize=(plot_ncol * sub_axes_unit, plot_nrow * sub_axes_unit))
        k=0
        ref_ghist, ref_gbins = cal_hist_density(all_bg)
        ref_x = (ref_gbins[1:] + ref_gbins[:-1]) / 2
        for i in range(plot_nrow):
            for j in range(plot_ncol):
                if plot_nrow==1:
                    ax=axes[j]
                else:
                    if plot_ncol==1:
                        ax=axes[i]
                    else:
                        ax = axes[i, j]
                if k<need_axes_num:
                    cell_v=sort_datas[:,k]
                    cell=cells[k]
                    if info is not None:
                        ax.set_title(f'{info[k]}')
                    ax.plot(ref_x, ref_ghist, '-', label='All')
                    ghist, gbins = cal_hist_density(cell_v)
                    ax.plot((gbins[1:] + gbins[:-1]) / 2, ghist, '-', label=cell)
                    ax.legend(loc='upper right')
                else:
                    plt.delaxes(ax)
                k+=1
        plt.tight_layout()
        make_dir(os.path.dirname(fig_path))
        plt.savefig(fig_path)

    def __plot_sc_depth(self,depth,n_gene,path):
        joint = seaborn.jointplot(x=depth, y=n_gene, kind="kde")
        joint.set_axis_labels('UMI count','Gene count')
        plt.title(f'Median: UM={np.median(depth)}, Gene={np.median(n_gene)}')
        plt.tight_layout()
        plt.savefig(path)

    def __plot_cell_umap(self,datas:[],cell_types,fig_path,labels=None):
        need_axes_num=len(datas)
        sub_axes_unit=3
        n = int(np.sqrt(need_axes_num) - 0.0001) + 1
        plot_ncol = n
        plot_nrow = int(need_axes_num/plot_ncol-0.0001)+1
        fig, axes = plt.subplots(plot_nrow, plot_ncol, figsize=(plot_ncol * sub_axes_unit, plot_nrow * sub_axes_unit))
        k=0
        for i in range(plot_nrow):
            for j in range(plot_ncol):
                if plot_nrow==1:
                    ax=axes[j]
                else:
                    ax = axes[i, j]
                if k<need_axes_num:
                    x,y=datas[k]
                    if labels is None:
                        ax.scatter(x, y, s=10)
                    else:
                        label_arr=np.array(labels[k])
                        unique_labels=np.unique(label_arr)
                        for label in unique_labels:
                            sidx=np.where(label_arr==label)
                            ax.scatter(x[sidx], y[sidx], label=label, s=10)
                    ax.set_title(f'{cell_types[k]}')
                    ax.legend()
                else:
                    plt.delaxes(ax)
                k+=1
        plt.tight_layout()
        make_dir(os.path.dirname(fig_path))
        plt.savefig(fig_path)

    def __plot_sample_size_k_scatter(self,cell_size,k_factors):
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.plot(cell_size, k_factors, '.')
        ax.set_xlabel('Sample size')
        ax.set_ylabel('k')
        plt.tight_layout()
        plt.savefig(f'{self.intermediate_dir}/cell_size-k-scatter.png')


    def load_expr(self,include_cells=None):
        is_expr_file=True
        if "expression_dir" in self.para and self.para["expression_dir"] is not None:
            is_expr_file=False
        if is_expr_file:
            exprs = self.__read_expr_file(include_cells)
        else:
            exprs=self.__read_expr_dir(include_cells)
        for expr in exprs:
            yield expr

    def load_cell_k(self):
        k_factor_path = f'{self.intermediate_dir}/k_factor.txt'
        if os.path.isfile(k_factor_path):
            cell_k = {}
            k_df=pd.read_table(k_factor_path,index_col=0,header=0)
            for c in k_df.columns:
                cell_k[c]=k_df.loc['k_factor',c]
            return cell_k
        else:
            return None

    def cal_expr_mean(self):
        out_path=f'{self.node_score_dir}/expr.mean.txt.gz'
        if self.para['keep_expr'] and os.path.isfile(out_path):
            log('The mean expression profile of genes has already been generated.')
            return
        cell_k=self.load_cell_k()
        if cell_k is None:
            log('No k factor file, can not generate mean expression.')
            return
        include_cells=sorted(cell_k.keys())
        expr_datas = self.load_expr(include_cells)
        mean_expr = []
        cells = []
        for expr, cell_name, gene_id in expr_datas:
            cells.append(cell_name)
            expr=self.__normalize_expr(expr)
            if self.para['log_trans_expr']:
                expr=np.log2(expr+1)
            mean_expr.append(np.nanmean(expr, axis=1))
        genes=gene_id['gene_symbol']
        df = pd.DataFrame(mean_expr, columns=genes, index=cells).T
        df.index.name = 'Gene'
        df.to_csv(out_path, sep='\t', float_format='%.4f')
        self.plot_expr_dist()
        log('Completed mean expression profile of genes.')

    def cell_type_qc(self,sum_degree:np.ndarray):
        data=sum_degree
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        threshold = 1.5 #1.5
        remain_index = np.where((data >= Q1 - threshold * IQR) & (data <= Q3 + threshold * IQR))[0]
        rmidx=[int(i) for i in remain_index]
        return rmidx

    def test_qc_cells(self):
        resample_path=f'{self.intermediate_dir}/resample_corr.txt.gz'
        df=pd.read_table(resample_path)
        cells=df.columns.values
        resample_corr=df.values
        ## filter cell types by global degree.
        mean_degree = np.nanmean(resample_corr, axis=0)
        print({cells[i]: mean_degree[i] for i in range(len(mean_degree))})
        remain_idx = self.cell_type_qc(mean_degree)
        remove_idx = np.array([i for i in range(len(cells)) if i not in remain_idx])
        log(f'filter out {len(remove_idx)} cell types: {cells[remove_idx]}, remain {len(remain_idx)}.')
        resample_corr = resample_corr[:, remain_idx]
        cells = cells[remain_idx]
        k_factors = self.__cal_k_factor(resample_corr)
        cell_k = {cells[i]: k_factors[i] for i in range(len(k_factors))}
        print(cell_k)
        log(f'calculated k factors')

    def cal_degree_adjusted(self):
        min_k=self.para['min_k']
        max_k=self.para['max_k']
        net_cate=['raw','k_adjusted'] # Notice: the order can not be changed.
        # net_cate=['k_adjusted']
        k_factor_path = f'{self.intermediate_dir}/k_factor.txt'
        degree_path = f'{self.node_score_dir}/degree.k_adjusted.txt.gz'
        resample_path=f'{self.intermediate_dir}/resample_corr.txt.gz'
        if self.para['keep_expr'] and os.path.isfile(k_factor_path):
            cell_k=self.load_cell_k()
            log('The k factor has already been calculated.')
        else:
            resample_size=self.para['resample_size']
            expr_datas = self.load_expr() #['Ast2','ExN9','Oli3']
            sample_size = []
            cells = []
            resample_corr=[]
            for expr, cell_name, gene_id in expr_datas:
                sample_size.append(expr.shape[1])
                seq_depth=np.nansum(expr,axis=0)
                cor_mat=np.abs(self.corr_matrix(expr,seq_depth))
                triu_mat = cor_mat[np.triu_indices(cor_mat.shape[0], k=1)]
                if len(triu_mat)<resample_size:
                    resample_size=len(triu_mat)
                np.random.seed(10086)
                random_resample = np.random.choice(triu_mat[~np.isnan(triu_mat)], resample_size, replace=False)
                resample_corr.append(random_resample)
                cells.append(cell_name)
                log(f'resampled from {cell_name}')
            sample_size = np.array(sample_size)
            cells=np.array(cells)
            raw_cell_num=len(cells)
            raw_cells=np.copy(cells)
            resample_corr=np.array(resample_corr).T
            rc_df=pd.DataFrame(resample_corr,columns=cells)
            rc_df.to_csv(resample_path,sep='\t',float_format='%.4f',index=None)
            ## filter cell types by global degree.
            mean_abs_r=np.nanmean(resample_corr,axis=0)
            remain_idx=self.cell_type_qc(mean_abs_r)
            mean_abs_r=mean_abs_r[remain_idx]
            resample_corr=resample_corr[:, remain_idx]
            sample_size=sample_size[remain_idx]
            cells=cells[remain_idx]
            tmp_resample_corr=np.copy(resample_corr)
            sec_remain_idx=np.arange(0,len(remain_idx))
            k_factors={}
            while True:
                continue_do=False
                tmp_k_factors = self.__cal_k_factor(tmp_resample_corr[:,sec_remain_idx])
                tmp_idx=[]
                for pi in range(len(tmp_k_factors)):
                    kv=tmp_k_factors[pi]
                    if kv>min_k and kv<max_k:
                        tmp_idx.append(sec_remain_idx[pi])
                    else:
                        continue_do=True
                sec_remain_idx=np.array(tmp_idx)
                if not continue_do:
                    k_factors=tmp_k_factors
                    break
            mean_abs_r = mean_abs_r[sec_remain_idx]
            sample_size = sample_size[sec_remain_idx]
            cells = cells[sec_remain_idx]
            remove_num=raw_cell_num-len(sec_remain_idx)
            remove_cells=[c for c in raw_cells if c not in cells]
            if remove_num>0:
                log(f'filter out {remove_num} cell types: {remove_cells}, remain {len(sec_remain_idx)}.')
            log(f'calculated k factors')
            cell_k={cells[i]:k_factors[i] for i in range(len(k_factors))}
            self.cell_k=cell_k
            ck_df=pd.DataFrame({'k_factor':k_factors,'sample_size':sample_size,'mean_abs_r':mean_abs_r},index=cells).T
            ck_df.index.name='index'
            ck_df.to_csv(k_factor_path,sep='\t')
            # info = [f'sample size={ss}' for ss in sample_size]
            # for cate in net_cate:
                # if cate=='k_adjusted':
                #     resample_corr=np.power(resample_corr,k_factors)
                # self.__plot_distribution(resample_corr,cells,f'{self.intermediate_dir}/edge_weight.dist.{cate}.png',info)
            self.__plot_sample_size_k_scatter(sample_size,k_factors)
            log(f'saved and plot k factors')
        if self.para['keep_expr'] and os.path.isfile(degree_path):
            log('The degree of gene co-expression network has already been calculated.')
        else:
            cells=[]
            sample_size = []
            cate_cell_degree={cate:[] for cate in net_cate}
            # only cells in k factor file
            include_cells=sorted(cell_k.keys())
            expr_datas = self.load_expr(include_cells)
            cell_depth_gene=[]
            for expr, cell_name, gene_id in expr_datas:
                seq_depth=np.sum(expr,axis=0)
                n_gene=np.count_nonzero(expr,axis=0)
                for s in range(len(seq_depth)):
                    cell_depth_gene.append([seq_depth[s],n_gene[s]])
                cells.append(cell_name)
                sample_size.append(expr.shape[1])
                cor_mat=np.abs(self.corr_matrix(expr,seq_depth))
                cor_mat[np.diag_indices(cor_mat.shape[0])] = 0
                for cate in net_cate:
                    if cate=='k_adjusted':
                        cor_mat=np.power(cor_mat,cell_k[cell_name])
                    cate_cell_degree[cate].append(np.nansum(cor_mat,axis=1))
                log(f'calculated degree of {cell_name}')
            genes = gene_id['gene_symbol']
            info = [f'sample size={ss}' for ss in sample_size]
            for cate in cate_cell_degree.keys():
                deg_vals=np.array(cate_cell_degree[cate]).T
                deg_df=pd.DataFrame(deg_vals,columns=cells,index=genes)
                deg_df.index.name='Gene'
                deg_df.to_csv(f'{self.node_score_dir}/degree.{cate}.txt.gz',sep='\t',float_format='%.4f')
                self.__plot_distribution(deg_vals,cells,f'{self.intermediate_dir}/degree.dist.{cate}.png',info)
            cdg=np.array(cell_depth_gene)
            self.__plot_sc_depth(cdg[:,0],cdg[:,1],f'{self.intermediate_dir}/sc_depth.png')
            log(f'saved and plot degree')

    def __read_df_sort_col(self,path, **argv):
        df = pd.read_table(path, **argv)
        cols = df.columns
        suffixs = set()
        cells = set()
        for c in cols:
            arr = str(c).split('.')
            if len(arr) > 1:
                suffixs.add(arr[-1])
                cells.add('.'.join(arr[:-1]))
            if len(arr) == 1:
                cells.add(c)
        cells = sorted(cells)
        if 'z' in suffixs:
            df = df[[f'{c}.z' for c in cells]]
            df.columns = cells
        else:
            df = df[cells]
        return df
    def expr_analysis(self):
        '''
        Analysis the relationship between expression and network.
        :return:
        '''
        kggsee_jar=f'{self.para["kggsee_dir"]}/kggsee.jar'
        kggsee_resources=f'{self.para["kggsee_dir"]}/resources'
        ## cal rez and plot corr of deg (and expr)
        expr_keys=['degree.k_adjusted']
        if self.para['run_expr_dese']:
            expr_keys.append('expr.mean')
        expr_dfs={}
        for ek in expr_keys:
            cmd=kggsee_rez(f'{self.node_score_dir}/{ek}.txt.gz',f'{self.intermediate_dir}/{ek}',3,kggsee_jar,kggsee_resources)
            run_command(cmd,'\t')
            df=self.__read_df_sort_col(f'{self.intermediate_dir}/{ek}.REZ.webapp.genes.txt.gz',header=0,index_col=0)
            df=df.loc[~df.index.duplicated(),:]
            expr_dfs[ek]=df
            corr=df.corr(method='spearman')
            corr=cluster_df(corr)
            cell_corr_heatmap(corr,f'{self.intermediate_dir}/{ek}.REZ.corr.png')

        ## compare deg and expr if need expr.
        if self.para['run_expr_dese']:
            ex_df=expr_dfs['expr.mean']
            de_df=expr_dfs['degree.k_adjusted']
            common_celltypes=sorted(set(ex_df.columns).intersection(de_df.columns))
            common_genes=sorted(set(ex_df.index).intersection(de_df.index))
            ex_df=ex_df.loc[common_genes,:]
            de_df=de_df.loc[common_genes,:]
            corr_coefs=[]
            for cc in common_celltypes:
                stat,p=scipy.stats.spearmanr(ex_df.loc[:,cc],de_df.loc[:,cc])
                corr_coefs.append([stat,p])
            rdf=pd.DataFrame(corr_coefs,index=common_celltypes,columns=['spearman_r','p'])
            rdf.to_csv(f'{self.intermediate_dir}/celltype.cross.REZ.corr.txt',sep='\t')
        pass

    def plot_expr_dist(self):
        expr_path=f'{self.node_score_dir}/expr.mean.txt.gz'
        df=pd.read_table(expr_path, header=0, index_col=0)
        expr_vals=df.values
        cells=df.columns.values.tolist()
        count = np.sum(expr_vals > 0.1, axis=0)
        data={}
        for i in range(len(count)):
            data[cells[i]]=count[i]
        for c in sorted(data,key=lambda x:-data[x]):
            print(f'{c}: {data[c]} ({expr_vals.shape[0]})')
        self.__plot_distribution(expr_vals, cells, f'{self.intermediate_dir}/expr.dist.png')

    def gene_based_analysis(self):
        '''
        run DESE framework using gene degree by KGGSEE toolkit.
        :return:
        '''
        score_files = []
        for f in os.listdir(self.node_score_dir):
            if f.startswith('degree.raw'):
                continue
            if not self.para['run_expr_dese'] and f.startswith('expr.mean'):
                continue
            score_files.append(f)
        score_files = sorted(score_files)
        gene_score_files=' '.join([f'"{self.node_score_dir}/{sf}"' for sf in score_files])
        common_parars=[
            self.para['gwas_path'],
            gene_score_files,
            f'{self.result_dir}/{self.para["phenotype_abbr"]}',
            f'{self.para["kggsee_dir"]}/kggsee.jar',
            f'{self.para["kggsee_dir"]}/resources',
            self.para["multiple_testing"],
            self.para["p_value_cutoff"],
            self.para["top_n_gene"],
            self.para["nt"],
            self.para["chrom_col"],
            self.para["pos_col"],
            self.para["p_col"],
            self.para["buildver"],
            self.para["rm_hla"],
            self.para["java_path"],
            self.para["jvm_gb"],
        ]
        if "vcf_ref" in self.para and self.para["vcf_ref"] is not None:
            pop=os.path.basename(self.para["vcf_ref"]).split('.')[0]
            # cmd=kggsee_dese(*common_parars,vcf_ref=self.para["vcf_ref"],
            #                 keep_ref=f"{self.saved_ref_dir}/{self.para['buildver']}/{pop}")
            cmd = kggsee_dese(*common_parars, vcf_ref=self.para["vcf_ref"])
        else:
            cmd=kggsee_dese(*common_parars,saved_ref=self.para["saved_ref"])
        log(f'{"*"*5} start to run gene based analysis function of KGGSEE toolkit. {"*"*5}')
        print(f'Command: {cmd}')
        print(f'KGGSEE Log:')
        code=run_command(cmd,'\t')
        log(f'{"*"*5} complete to run KGGSEE. {"*"*5}')
        return code

    def load_k_factor(self):
        cell_k = {}
        k_factor_path = f'{self.intermediate_dir}/k_factor.txt'
        k_df = pd.read_table(k_factor_path, index_col=0, header=0)
        for c in k_df.columns:
            cell_k[c] = k_df.loc['k_factor', c]
        return cell_k

    def detect_gene_module(self):
        '''
        Detect disease associated gene network modules.
        :return:
        '''
        n_top_genes = 1000000
        gene_score_n_top_genes = self.para['module_gene_score_n_top_genes']
        expr_norm_method = self.para['normalize_expr']
        cut_abs_r = self.para['module_cut_edge_weight']
        assoc_cell_p_cutoff=self.para['assoc_cell_p_cutoff']
        gene_score_file=None
        if self.para['run_expr_dese']:
            gene_score_file=self.node_score_file_name
        kggsee_prefix=f'{self.result_dir}/{self.para["phenotype_abbr"]}'
        kg = KGGSEE(kggsee_prefix)
        node_top_k = gene_score_n_top_genes
        # node_top_perc = 0.25
        module_node_size = 300
        module_min_gene = 20
        min_top_cell_num = 1 ## min number of cell to detect gene modules.

        cell_k=self.load_k_factor()
        assoc_cells=kg.assoc_sig_cells(assoc_cell_p_cutoff,gene_score_file,min_top_cell_num)
        assoc_genes=kg.cond_sig_assoc_gene(gene_score_file)
        cell_gene_score = {}
        module_net = {}
        if len(assoc_cells)==0:
            log(f'no significant cell types, skip detecting modules.')
            return None
        log(f'start to detect associated modules in {len(assoc_cells)} cell types: {", ".join(assoc_cells)}.')
        expr_datas=self.load_expr(assoc_cells)
        dfs=[]
        for expr, cell_name, gene_id in expr_datas:
            seq_depth=np.nansum(expr,axis=0)
            ## select only top expressed genes
            include_genes = None
            print(expr,gene_id,n_top_genes,expr_norm_method,include_genes)
            expr, gene_id = self.extract_expr_top(expr,gene_id,n_top_genes,expr_norm_method,include_genes)
            log(f'{cell_name}: load data and cal adj matrix')
            genes=gene_id['gene_symbol']
            raw_cor_mat=self.corr_matrix(expr,seq_depth)
            cor_mat=np.abs(raw_cor_mat)
            sign_cor_mat=raw_cor_mat/cor_mat
            if self.node_score_file_name=='degree.k_adjusted.txt.gz':
                cor_mat=np.power(cor_mat,cell_k[cell_name])
            cor_mat[np.diag_indices(cor_mat.shape[0])] = 0
            degree=np.nansum(cor_mat,axis=0)
            non_zero_idx=np.where(degree>0)[0]
            genes=genes[non_zero_idx]
            cor_mat=cor_mat[non_zero_idx,:][:,non_zero_idx]
            sign_cor_mat=sign_cor_mat[non_zero_idx,:][:,non_zero_idx]
            assoc_gene_idxs=np.array([i for i in range(len(genes)) if genes[i] in assoc_genes])
            ## gene score based sim_mat or cor_mat
            scores = np.nansum(cor_mat[assoc_gene_idxs, :], axis=0)
            cell_gene_score[cell_name] = scores
            log(f'{cell_name}: filter node')
            max_node_top_k = len(scores) - len(assoc_gene_idxs)
            if node_top_k > max_node_top_k:
                log(f'{cell_name}: set node_top_k={max_node_top_k}')
                node_top_k = max_node_top_k
            k = node_top_k
            copy_scores = np.copy(scores)
            copy_scores[assoc_gene_idxs] = -np.inf
            top_node_idx = np.argpartition(copy_scores, -k)[-k:]
            remain_node_idx = np.concatenate((top_node_idx, assoc_gene_idxs))
            ## exclude low connected nodes
            final_cor_mat=cor_mat[remain_node_idx, :][:, remain_node_idx]
            final_cor_mat[np.diag_indices(final_cor_mat.shape[0])] = 1
            final_sign_cor_mat=sign_cor_mat[remain_node_idx, :][:, remain_node_idx]
            log(f'{cell_name}: detect module')
            gn = GeneNetwork()
            modules = gn.recursive_detect_module(final_cor_mat, module_node_size, False,cut_abs_r)
            modules = [m for m in modules if len(m) >= module_min_gene]
            log(f'{cell_name}: test module')
            cell_minfo = []
            for i in range(len(modules)):
                m = modules[i]
                m_idxs = remain_node_idx[m]
                assoc_gene_in_m_idxs = sorted(set(m_idxs).intersection(assoc_gene_idxs))
                enrich_pv = stats.hypergeom.sf(len(assoc_gene_in_m_idxs), len(remain_node_idx), len(assoc_gene_idxs),
                                               len(m))
                module_id = f'{cell_name}_{i + 1}'
                cell_minfo.append([module_id, cell_name, enrich_pv, len(assoc_gene_in_m_idxs), len(m),
                                   ','.join(sorted(genes[assoc_gene_in_m_idxs])),
                                   ','.join(sorted(genes[m_idxs]))])
                cor_df=pd.DataFrame(final_cor_mat[m, :][:, m]*final_sign_cor_mat[m,:][:,m],columns=genes[m_idxs],
                                    index=genes[m_idxs])
                module_net[module_id] = {'cor_mat_df': cor_df}

            df = pd.DataFrame(cell_minfo, columns=['module_id', 'cell_type', 'p', 'assoc_gene_number',
                                             'module_gene_number',
                                             'assoc_gene', 'module_gene'])
            for pcol in ['p']:
                df[f'{pcol}.adj_fdr'] = pvalue_adjust(df[pcol].values.tolist())
            dfs.append(df)
        fdf=pd.concat(dfs,axis=0,ignore_index=True)
        fdf = fdf[['module_id','cell_type', 'p', 'p.adj_fdr','assoc_gene_number',
                 'module_gene_number', 'assoc_gene', 'module_gene']]
        fdf = fdf.sort_values(by=['cell_type','p'])
        fdf.to_csv(f'{kggsee_prefix}.assoc_gene_module.txt', sep='\t', index=False)
        with open(f'{kggsee_prefix}.assoc_gene_module.network.pyd', 'wb') as bw:
            pickle.dump(module_net, bw)

    def annotate_module(self):
        log(f'start annotate gene modules')
        kggsee_prefix=f'{self.result_dir}/{self.para["phenotype_abbr"]}'
        assoc_module_path = f'{kggsee_prefix}.assoc_gene_module.txt'
        function_enrich_p_cutoff = self.para['function_enrich_p_cutoff']
        function_enrich_top_n = self.para['function_enrich_top_n']
        anno_dbs=['GO:BP','GO:CC','GO:MF','KEGG']
        if not os.path.isfile(assoc_module_path):
            log(f'no assoc module file: {assoc_module_path}')
            return
        df = pd.read_table(assoc_module_path, header=0,index_col=0)
        df['p.global_adj_fdr']=pvalue_adjust(df['p'])
        df = df.loc[df['p.adj_fdr'] <= self.module_p_cutoff, :]
        df = df.sort_values(by=['p'])
        mid_genes={}
        for mid in df.index:
            module_genes = df.loc[mid, 'module_gene'].split(',')
            mid_genes[mid]=module_genes
        gp = GProfiler(return_dataframe=True)
        anno_df=gp.profile(organism='hsapiens',query=mid_genes,sources=anno_dbs,user_threshold=function_enrich_p_cutoff)
        for mid in df.index:
            for db in anno_dbs:
                anno_text='.'
                sdf=anno_df.loc[(anno_df['query']==mid) & (anno_df['source']==db),:]
                if sdf.shape[0]>0:
                    sdf=sdf.sort_values(by=['p_value'])
                    anno_terms=[f"{sdf.loc[i, 'name']} ({sdf.loc[i, 'native']}, P={sdf.loc[i,'p_value']:.2g})" for i in sdf.index]
                    show_n=len(anno_terms)
                    if show_n>function_enrich_top_n:
                        show_n=function_enrich_top_n
                    anno_text=';'.join(anno_terms[:show_n])
                df.loc[mid,db]=anno_text
        df.index.name = 'module_id'
        assoc_module_anno_path = f'{kggsee_prefix}.assoc_gene_module_anno.xlsx'
        df.to_excel(assoc_module_anno_path)

    def __get_gene_pvalue(self,genes:[]):
        ## associated gene p values. return [p,...], non-sig gene with p=1.
        gene_score_file=None
        if self.para['run_expr_dese']:
            gene_score_file=self.node_score_file_name
        kg_prefix=f'{self.result_dir}/{self.para["phenotype_abbr"]}'
        kg=KGGSEE(kg_prefix)
        gene_p=kg.cond_sig_assoc_gene_p(gene_score_file)
        pvs=[]
        for g in genes:
            p=1
            if g in gene_p:
                p=gene_p[g]
            pvs.append(p)
        return pvs

    def plot_module(self,include_mid=None):
        module_plot_cut_abs_r = self.para['module_plot_cut_edge_weight']
        top_n_degree = self.para['show_node_label_top_n']
        log(f'start plot gene modules')
        kggsee_prefix=f'{self.result_dir}/{self.para["phenotype_abbr"]}'
        assoc_module_path = f'{kggsee_prefix}.assoc_gene_module.txt'
        module_network_path=f'{kggsee_prefix}.assoc_gene_module.network.pyd'
        module_plot_dir=f'{kggsee_prefix}.plot'
        make_dir(module_plot_dir)
        if not os.path.isfile(assoc_module_path) or not os.path.isfile(module_network_path):
            return
        module_cor_mat={}
        with open(module_network_path, 'rb') as br:
            module_cor_mat=pickle.load(br)
        # cmap_heatmap=heatmap_color_rstyle_single()
        cmap_heatmap = heatmap_color_rstyle()
        cmap_heatmap_single = heatmap_color_rstyle_single()
        df = pd.read_table(assoc_module_path, header=0,index_col=0)
        df['p.global_adj_fdr']=pvalue_adjust(df['p'])
        df = df.loc[df['p.adj_fdr'] <= self.module_p_cutoff, :]
        df = df.sort_values(by=['p'])
        mid_genes={}
        for mid in df.index:
            if include_mid is not None and mid not in include_mid:
                continue
            cor_df=module_cor_mat[mid]['cor_mat_df']
            m_genes=cor_df.index.values.tolist()
            mid_genes[mid]=m_genes
            p_values = self.__get_gene_pvalue(m_genes)
            sig_genes=[]
            for i in range(len(m_genes)):
                if p_values[i]<1:
                    sig_genes.append(m_genes[i])
            raw_cor_mat=cor_df.values
            cor_mat=np.abs(raw_cor_mat)
            cor_mat[np.diag_indices(cor_mat.shape[0])] = 0
            degree = np.nansum(cor_mat, axis=0)
            node_size = np.power(5, degree / np.max(degree)) * 2 + 0.5
            node_size=node_size*3
            nodes = []
            links = []
            categories = [opts.GraphCategory(name='UnSig. genes'),
                          opts.GraphCategory(name='Sig. genes')]
            if top_n_degree>len(degree):
                top_n_degree=len(degree)
            top_degree_cut=sorted(degree)[-top_n_degree]
            for i in range(cor_mat.shape[0]):
                if not np.any(cor_mat[i, :] >= module_plot_cut_abs_r):
                    continue
                label_color='black'
                font_weight='normal'
                category = 0
                show = True
                # if degree[i] >= per_cut or m_genes[i] in sig_genes:
                #     show = True
                if degree[i] >= top_degree_cut:
                    label_color='red'# '#FC8452'
                    font_weight='bolder'
                if p_values[i] < 1:
                    category = 1
                nodes.append(opts.GraphNode(name=m_genes[i],symbol_size=node_size[i],category=category,value=degree[i],
                                            # itemstyle_opts=opts.ItemStyleOpts(color=node_color),
                                            label_opts=opts.LabelOpts(is_show=show,color=label_color,font_weight=font_weight)))
            nodes = sorted(nodes, key=lambda x: x.opts['category'])
            for i in range(cor_mat.shape[0] - 1):
                for j in range(i + 1, cor_mat.shape[0]):
                    w = cor_mat[i, j]
                    if w < module_plot_cut_abs_r:
                        continue
                    edge_color='#E1E1E1'
                    if m_genes[i] in sig_genes or m_genes[j] in sig_genes:
                        edge_color='#EE6666'
                    # links.append({"source": m_genes[i], "target": m_genes[j], "value": w})
                    links.append(opts.GraphLink(source=m_genes[i],target=m_genes[j],value=w,
                                                linestyle_opts=opts.LineStyleOpts(color=edge_color,curve=0.2,width=0.7)))
            graph = (
                Graph(init_opts=opts.InitOpts(width="1200px", height="800px"))
                .add("", nodes, links, categories, layout='circular',  # circular,force
                     is_draggable=True,
                     repulsion=10,
                     # linestyle_opts=opts.LineStyleOpts(width=0.5, curve=0.2),  # 设置边的样式
                     # label_opts=opts.LabelOpts(font_family='Arial', color='black'),
                     is_rotate_label=True,
                     )
                .set_global_opts(title_opts=opts.TitleOpts(title=f"{mid}"),legend_opts=opts.LegendOpts(
                    textstyle_opts=opts.TextStyleOpts(font_family='Arial',color='black'),orient='vertical'))
                .set_colors(['#73C0DE','#EE6666'])
            )

            graph.render(f"{module_plot_dir}/{mid}.module_network.html")
            cor_mat[np.diag_indices(cor_mat.shape[0])] = 1
            # fsize=10/40*cor_df.shape[0]
            seaborn.clustermap(cor_df, cmap=cmap_heatmap,vmin=-1,vmax=1)
            plt.savefig(f"{module_plot_dir}/{mid}.module_heatmap_sign.png")
            plt.close()
            seaborn.clustermap(cor_df.abs(), cmap=cmap_heatmap_single,vmin=0,vmax=1)
            plt.savefig(f"{module_plot_dir}/{mid}.module_heatmap.png")
            plt.close()

        # plot module similarity
        mids=sorted(mid_genes.keys())
        sim_mat=np.zeros(shape=[len(mids),len(mids)])
        for i in range(len(mids)):
            for j in range(i,len(mids)):
                ji=jaccard(mid_genes[mids[i]],mid_genes[mids[j]])
                sim_mat[i,j]=ji
                sim_mat[j,i]=ji
        sim_df=pd.DataFrame(sim_mat,columns=mids,index=mids)
        if sim_df.shape[0]>1:
            cmap=heatmap_color_rstyle_single()
            seaborn.clustermap(sim_df,cmap=cmap, annot=True, fmt=".2f",linewidths=0.5, linecolor='gray')
            plt.savefig(f"{module_plot_dir}/module_jaccard.png")
            plt.close()

    def plot_module_enrich_result(self,include_mid=None):
        log(f'start plot gene modules enrichment')
        kggsee_prefix=f'{self.result_dir}/{self.para["phenotype_abbr"]}'
        assoc_module_path = f'{kggsee_prefix}.assoc_gene_module.txt'
        module_plot_dir=f'{kggsee_prefix}.plot'
        make_dir(module_plot_dir)
        if not os.path.isfile(assoc_module_path):
            return
        df = pd.read_table(assoc_module_path, header=0,index_col=0)
        df['p.global_adj_fdr']=pvalue_adjust(df['p'])
        df = df.loc[df['p.adj_fdr'] <= 0.1, :]
        df = df.sort_values(by=['p'])
        for mid in df.index:
            if include_mid is not None and mid not in include_mid:
                continue
            mgenes=str(df.loc[mid,'module_gene']).strip().split(',')
            out_fig_path=f'{module_plot_dir}/{mid}.module_go_enrich.png'
            go_enrich_plot(mgenes,out_fig_path)


    def main(self):
        log('****** DGN ******')
        self.__print_effect_para()
        log('start ')
        self.cal_degree_adjusted()
        if self.para['run_expr_dese']:
            self.cal_expr_mean()
        # self.expr_analysis()
        self.gene_based_analysis()
        self.detect_gene_module()
        self.annotate_module()
        self.plot_module()
        log('finish ')

    def main_expr(self):
        log('****** DGN ******')
        self.__print_effect_para()
        log('start process expr and network')
        self.cal_degree_adjusted()
        if self.para['run_expr_dese']:
            self.cal_expr_mean()
        self.expr_analysis()
        log('finish degree step')

    def main_gwas(self):
        log('****** DGN ******')
        self.__print_effect_para()
        log('start gene and module based analysis')
        self.gene_based_analysis()
        log('finish gene-based step')

    def main_module(self):
        log('****** DGN ******')
        self.__print_effect_para()
        log('start module based analysis')
        self.detect_gene_module()
        self.annotate_module()
        self.plot_module()
        log('finish module step')


def trans_options(args):
    para=vars(args)
    return para

def main():
    parser = argparse.ArgumentParser(description='DGN: a Disease-associated Gene Network module analysis framework.')
    gene_expr=parser.add_argument_group('Input 1: Gene expression profile','Input gene expression matrix from single cell RNA-seq.')
    gene_expr.add_argument('--expression_dir',type=str,required=True,help='Input a directory that contains gene expression profile of different cell types.'
                                                                  'Each file represents a specific cell type, with filenames following the convention "{CellTypeName}.txt.gz".')

    gwas=parser.add_argument_group('Input 2: GWAS summary statistics','')
    gwas.add_argument('--gwas_path',type=str,help='GWAS summary table with header.',required=True)
    gwas.add_argument('--chrom_col', type=str, help='Chromosome column in GWAS summary table.', default='CHR')
    gwas.add_argument('--pos_col', type=str, help='Base position column in GWAS summary table.', default='BP')
    gwas.add_argument('--p_col', type=str, help='P-value column in GWAS summary table.', default='P')
    gwas.add_argument('--buildver', type=str,default='hg19',choices=['hg19','hg38'],
                      help='Specifies the reference genome version of the coordinates.')
    gwas.add_argument('--phenotype_abbr', type=str, help='Phenotype abbreviation, used as the file identifier for this phenotype in the output results', required=True)

    output=parser.add_argument_group('Output', '')
    output.add_argument('--output_dir',type=str,help='Output directory.',required=True)

    step_1=parser.add_argument_group('Step 1: Construction of gene co-expression network', '')
    step_1.add_argument('--trans_gene_symbol',action='store_true',help='Convert the gene ID in the expression profiles to gene symbol of HGNC. ')
    step_1.add_argument('--min_expr_value', type=float, default=0,
                           help='The minimum average gene expression level used for constructing co-expression networks.')
    step_1.add_argument('--min_cell', type=int, default=100,help='The minimum number of cells or samples required in the analysis.')
    step_1.add_argument('--max_cell', type=int, default=1000,
                           help='The max number of cells or samples required in the analysis.')
    step_1.add_argument('--min_k', type=float, default=0.5,help='The minimum value of k in normalization for co-expression network.')
    step_1.add_argument('--max_k', type=float, default=1.5,help='The max value of k in normalization for co-expression network.')
    step_1.add_argument('--edge_method', type=str, default='cs-core',choices=['pearson','cs-core'],
                           help='The method for calculating edge weights (i.e., gene correlations)')
    step_1.add_argument('--resample_size', type=int, default=100000,
                           help='Sample size for resampling edge weights when correcting the co-expression network.')
    step_1.add_argument('--keep_expr', action='store_true',
                           help='No need to recalculate the gene centrality in the gene co-expression network for the next analysis.')

    step_2=parser.add_argument_group('Step 2: Disease-associated genes and cell types', '')
    # step_2.add_argument('--kggsee_dir',type=str,help='Directory of downloaded and extracted KGGSEE.',required=True)
    step_2.add_argument('--vcf_ref', type=str,required=True, help='Specifies a VCF file of genotypes sampled from a reference population.')
    step_2.add_argument('--multiple_testing', type=str,default='benfdr',choices=['benfdr','bonf','fixed'],
                        help='Specifies the method for multiple testing correction.')
    step_2.add_argument('--p_value_cutoff', type=float, default=0.05,
                        help='Specifies the threshold of the adjusted p-value for fine-mapping.')
    step_2.add_argument('--top_n_gene', type=int, default=1000,
                        help='Maximum number of genes entering conditional association analysis.')
    step_2.add_argument('--nt', type=int, default=1,
                        help='Specifies the number of threads.')
    step_2.add_argument('--rm_hla', action='store_true',help='Remove HLA region.')
    step_2.add_argument('--java_path', type=str, default='java', help='Java path.')
    step_2.add_argument('--jvm_gb', type=int, default='40', help='JVM value (GB).')

    step_3=parser.add_argument_group('Step 3: Disease-associated gene network modules', '')
    step_3.add_argument('--assoc_cell_p_cutoff', default=0.05, type=float,
                           help='The adjusted p cutoff for associating cell types. The significant cell types will be used for associated network module analysis.')
    step_3.add_argument('--module_gene_score_n_top_genes', type=int, default=5000,
                           help='In the module detection function, specify the number of top genes with high disease-related scores for module detection')
    step_3.add_argument('--module_cut_edge_weight', type=float, default=0.5,
                           help='Minimum edge weight for pruning nodes in modules.')
    step_3.add_argument('--module_plot_cut_edge_weight', type=float, default=0.6,
                           help='Minimum edge weight for pruning nodes in plotting modules.')
    step_3.add_argument('--show_node_label_top_n', type=int, default=5,
                           help='Display the labels of the top n genes by connectivity.')
    step_3.add_argument('--function_enrich_p_cutoff', type=float, default=0.05,
                           help='P-value threshold for functional enrichment analysis by g:Profiler')
    step_3.add_argument('--function_enrich_top_n', type=int, default=3,
                           help='Maximum number of functional enrichment terms displayed.')

    step_extra=parser.add_argument_group('Extra analysis: DESE, inferring disease-associated genes and cell types using '
                                     'gene expression.', '')
    step_extra.add_argument('--run_expr_dese', action='store_true',
                           help='Calculate the mean expression profile of genes and run DESE.')
    step_extra.add_argument('--normalize_expr', type=str, default='no',choices=['no','cpm'],
                           help='The normalization method for the expression profiles')
    step_extra.add_argument('--log_trans_expr', action='store_true',
                           help='Transform the expression value into log2(x+1).')

    func=parser.add_argument_group('Separated functions', '')
    func.add_argument('--detect_module', action='store_true',
                           help='Detect module function.')

    parser.parse_args()
    args = parser.parse_args()
    para=trans_options(args)

    dgn=DGN(para)
    if para['detect_module']:
        dgn.main_module()
    else:
        dgn.main()

if __name__ == '__main__':
    if len(sys.argv)==2 and sys.argv[1]=='--version':
        print(f'Congratulations! {PKG_NAME} {VERSION} is installed successfully. \nRun `python dgn.py -h` to view usage instructions.')
    else:
        main()

