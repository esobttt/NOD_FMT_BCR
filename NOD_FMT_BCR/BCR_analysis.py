import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
import numpy as np
import math
import scipy.stats
from matplotlib.ticker import MaxNLocator
import networkx as nx
from itertools import combinations
from matplotlib.transforms import TransformedBbox
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
from matplotlib import colors
from matplotlib.ticker import FormatStrFormatter
import cblind as cb
from scipy.stats import ttest_ind, shapiro, mannwhitneyu, levene, chi2_contingency, kendalltau, wilcoxon

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 20

color_list = ['#2166AC', '#B2182B', '#5D3A9B']
color_list = ['#B2182B', '#2166AC', '#5D3A9B']
font = 'Arial'
figsize = (7, 5)

savedir = '/home/ajhml71531/sjogren_mouse/figure/Figure4'

isotypes = ['M', 'D', 'A', 'G1', 'G2B', 'G2C', 'G3', 'E']

target_organ = 'LG'
target_exp = 'FMT1'

if target_exp == 'FMT1':
    samples = ['NOD-FMT', 'PBS-T']
else:
    samples = ['HC-FMT', 'SS-FMT']
    
sample_color = {samples[i]: color_list[i] for i in range(len(samples))}
    
    
def select_significance(list_1,list_2):
    if (len(list_1) < 5) or (len(list_2) < 5):
        return

    if (shapiro(list_1)[1] > 0.05) and (shapiro(list_2)[1] > 0.05):
        if levene(list_1, list_2)[1] > 0.05:
            t,p = ttest_ind(list_1,list_2, equal_var=True)
            return 't-test_ind', p
        else:
            t,p = ttest_ind(list_1,list_2, equal_var=False)
            return 't-test_welch', p
    
    else:
        t,p = mannwhitneyu(list_1,list_2)
        return 'Mann-Whitney', p

## group-specific shared clonotype comparison
if True:
    if target_organ == 'LG':
        shared_dir = '/home/ajhml71531/sjogren_mouse/FMT/contam_elimination/analysis/'
    else:
        shared_dir = '/home/ajhml71531/sjogren_mouse/FMT/contam_elimination_SP/analysis/'
    if target_exp == 'FMT1':
        shared_file_1 = os.path.join(shared_dir, 'FMT1_E_shared_exclude_FMT1_P_FMT1_N.tsv')
        shared_file_2 = os.path.join(shared_dir, 'FMT1_P_shared_exclude_FMT1_E_FMT1_N.tsv')
    else:
        shared_file_1 = os.path.join(shared_dir, 'FMT2_E_shared_exclude_FMT2_P_FMT_N.tsv')
        shared_file_2 = os.path.join(shared_dir, 'FMT2_P_shared_exclude_FMT2_E_FMT_N.tsv')

    shared_files = [shared_file_1, shared_file_2]

    conv_dict = {k: {} for k in samples}
    conv_dict2 = {'group': [], 'degree of sharing': [], 'count': []}
    iso_dict = {s: {k: 0 for k in isotypes} for s in samples}
    degree_dict = {s: [] for s in samples}
    freq_dict = {s: [] for s in samples}
    div_dict = {s: [] for s in samples}
    sample_dict = {s: [] for s in samples}

    for i in range(len(samples)):
        cnt = 0
        with open(shared_files[i], 'r') as handle:
            header = handle.readline().strip().split('\t')
            col_dict = {}
            for h in header:
                col_dict[h] = header.index(h)
            while True:
                line = handle.readline()
                if not line:
                    break
                d = line.split('\t')
                d[-1] = d[-1].strip()
                convergence = int(d[col_dict['convergence']])
                max_freq = float(d[col_dict['freq_max']])
                div = float(d[col_dict['divergence']])
                num_seq = int(d[col_dict['num_sequences']])

                if convergence < 2:
                    continue
                    
                sample = [s.split('=')[0][1:] for s in d[col_dict['samples']][1:-1].split('|')]
                LG_cnt = 0
                for s in sample:
                    if 'LG' in s:
                        LG_cnt += 1
                
                if convergence not in conv_dict[samples[i]]:
                    conv_dict[samples[i]][convergence] = 1
                else:
                    conv_dict[samples[i]][convergence] += 1
                    
                
                degree_dict[samples[i]].append(convergence)
                freq_dict[samples[i]].append(max_freq)
                div_dict[samples[i]].append(div)
                sample_dict[samples[i]].append(sample)
                
                for iso in isotypes:
                    iso_dict[samples[i]][iso] += float(d[col_dict[iso]]) * num_seq
                
                cnt += 1
                
    total_freq = {}
    for s in samples:
        total_freq[s] = [np.log10(f) for f in freq_dict[s]]

    total_mean = [np.mean(freq_dict[s]) for s in samples]
    total_std = [np.std(freq_dict[s]) for s in samples]
    
    for j in range(len(samples)):
        if j > 0:
            continue
        for i in range(2, 11):
            if i not in conv_dict[samples[j]].keys():
                conv_dict[samples[j]][i] = 0
         
            
    for k in conv_dict:
        for d in conv_dict[k]:
            if (k == samples[1]) and (d == 10):
                continue
            conv_dict2['group'].append(k)
            conv_dict2['count'].append(conv_dict[k][d])
            if k == samples[0]:
                conv_dict2['degree of sharing'].append(d*10)
            else:
                conv_dict2['degree of sharing'].append(d*100/9)
            
    df = pd.DataFrame.from_dict(conv_dict2, orient='index').T
    
    
    
    ### samplewise heatmap    
    if True:
        w_df = pd.read_csv('/home/ajhml71531/sjogren_mouse/FMT/analysis/samplewise_sim.tsv', sep='\t')
        w_dict = {'from': w_df['from'].tolist() + w_df['to'].tolist(), 'to': w_df['to'].tolist() + w_df['from'].tolist(), 'sim': w_df['sim'].tolist()*2 }
        all_sample = list(set(w_dict['from']))
        for s in all_sample:
            w_dict['from'].append(s)
            w_dict['to'].append(s)
            w_dict['sim'].append(1)
        w_df = pd.DataFrame(w_dict)

        w_df = w_df[((w_df['from'].str.contains('FMT1')) | (w_df['to'].str.contains('FMT1')))]
        w_df = w_df[((~w_df['from'].str.contains('N')) & (~w_df['to'].str.contains('N')))]
        w_df = w_df.pivot('to', 'from', 'sim')
    
        w_df = w_df.reindex(['FMT1_E_LG%s' % str(i+1) for i in range(10)] + ['FMT1_P_LG%s' % str(i+1) for i in range(9)], axis=1)
        w_df = w_df.reindex(['FMT1_E_LG%s' % str(i+1) for i in range(10)] + ['FMT1_P_LG%s' % str(i+1) for i in range(9)], axis=0)
  
        groups = {'NOD-FMT': ['FMT1_E_LG%s' % str(i+1) for i in range(10)],
                  'PBS-T': ['FMT1_P_LG%s' % str(i+1) for i in range(9)],
                  }

        cmap = cb.cbmap('YlGnBu_r')
        fig, ax = plt.subplots(figsize=(12, 10))
        sns.heatmap(w_df, ax=ax, cmap=cmap, vmax=0.01)
        ax.hlines([10, 19], *ax.get_xlim(), colors='w', linewidths=3)
        ax.vlines([10, 19], *ax.get_ylim(), colors='w', linewidths=3)
        ax.set_ylabel("")
        ax.set_xlabel("")

        ax.set_xticklabels(['NOD-FMT%s' % str(i+1) for i in range(10)] + ['PBS-T%s' % str(i+1) for i in range(9)], fontsize=20)
        ax.set_yticklabels(['NOD-FMT%s' % str(i+1) for i in range(10)] + ['PBS-T%s' % str(i+1) for i in range(9)], fontsize=20)

        ax2 = ax.twiny()

        ax2.spines["bottom"].set_position(("axes", -0.25))
        ax2.tick_params('both', length=0, width=0, which='minor')
        ax2.tick_params('both', length=0, width=0,direction='in', which='major')
        ax2.xaxis.set_ticks_position("bottom")
        ax2.xaxis.set_label_position("bottom")

        ax2.set_xticks([0.0, 10/19, 1])
        ax2.xaxis.set_major_formatter(ticker.NullFormatter())
        ax2.xaxis.set_minor_locator(ticker.FixedLocator([10/38, 29/38]))
        ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['NOD-FMT', 'PBS-T']))

        ax3 = ax.twinx()

        ax3.spines["left"].set_position(("axes", -0.25))
        ax3.tick_params('both', length=0, width=0, which='minor')
        ax3.tick_params('both', length=0, width=0,direction='in', which='major')
        ax3.yaxis.set_ticks_position("left")
        ax3.yaxis.set_label_position("left")

        ax3.set_yticks([0.0, 9/19, 1])
        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
        ax3.yaxis.set_minor_locator(ticker.FixedLocator([9/38, 29/38]))
        ax3.yaxis.set_minor_formatter(ticker.FixedFormatter(['PBS-T', 'NOD-FMT']))

        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=20)
        plt.show()

    ## clonotype similarity boxplot 
    if True:
        w_df = pd.read_csv('/home/ajhml71531/sjogren_mouse/FMT/analysis/samplewise_sim.tsv', sep='\t')
        samples = ['FMT1_E', 'FMT1_P']
        plot_group_dict = {'FMT1_N': 'NC','FMT1_P': 'PBS-T', 'FMT1_E': 'NOD-FMT'}

        plot_group_name = [plot_group_dict[s] for s in samples]
        sample_sim_dict = {}
        target_idx = 0
        target = samples[target_idx]
        target = None

        target_pairs = [(samples[i], samples[i]) for i in range(len(samples))] + [(samples[0], samples[1])]
        for p in target_pairs:
         
            if p[0] == 'FMT_N':
                s_df = w_df[((w_df['from'].str.contains('N')) & (w_df['to'].str.contains(p[1]))) | ((w_df['to'].str.contains('N')) & (w_df['from'].str.contains(p[1])))]
            elif p[1] == 'FMT_N':
                s_df = w_df[((w_df['from'].str.contains('N')) & (w_df['to'].str.contains(p[0]))) | ((w_df['to'].str.contains('N')) & (w_df['from'].str.contains(p[0])))]
            else:
                s_df = w_df[((w_df['from'].str.contains(p[0])) & (w_df['to'].str.contains(p[1])) | (w_df['to'].str.contains(p[0])) & (w_df['from'].str.contains(p[1])))]
            if target == None:
    
                sample_sim_dict[p] = s_df['sim'].tolist()
            else:
                if target == p[0]:
                    sample_sim_dict[p[1]] = s_df['sim'].tolist()
                else:
                    sample_sim_dict[p[0]] = s_df['sim'].tolist()
        if target != None:
            for s in samples:
                if s != target:
                    continue
                if s == 'FMT_N':
                    s_df = w_df[((w_df['from'].str.contains('N')) & w_df['to'].str.contains('N'))]
                else:
                    s_df = w_df[((w_df['from'].str.contains(s)) & w_df['to'].str.contains(s))]
                sample_sim_dict[s] = s_df['sim'].tolist()
       
        sample_sim_df = pd.DataFrame.from_dict(sample_sim_dict, orient='index').T
    
        plot_sample_name = [plot_group_name[i] + ' vs ' + plot_group_name[i] for i in range(len(plot_group_name))] + [plot_group_name[0] + ' vs ' + plot_group_name[1]] 

        sample_sim_df.columns = plot_sample_name

        
        cb_color = cb.Colorplots().cblind(12)[0]
        cb_index = [10, 11]
        cb_index = [6, 4]
        cb_index = [10, 11, 6, 4]
       
        colors = {plot_sample_name [i]: cb_color[cb_index[i]] for i in range(len(plot_sample_name))}
        colors = {plot_sample_name [i]: color_list[i] for i in range(len(plot_sample_name))}
        pairs = [(plot_sample_name[i], plot_sample_name[target_idx]) for i in range(len(plot_group_name)) if i != target_idx]
        pairs = list(combinations(plot_sample_name , 2))
            
        fig, ax = plt.subplots(figsize=(14, 9))

        fig.text(0.03, 0.50, "Clonotype similarity", va='center', rotation = 'vertical')
        
        sns.boxplot(data=sample_sim_df, palette=colors, width=0.7, showfliers=False, boxprops=dict(alpha=.8))
        sns.stripplot(data=sample_sim_df, palette=colors, size=7)

        for i,artist in enumerate(ax.artists):
            col = color_list[i]
            artist.set_edgecolor(col)    

            for j in range(i*5,i*5+5):
                line = ax.lines[j]
                if (j %5 == 3) or (j%5 == 2):
                    line.set_color('white')
                else:
                    line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)
                if j % 5 == 4:
                    line.set_linewidth(3)
            
        for p in pairs:
            s, pp = select_significance(sample_sim_df[p[1]].dropna().tolist(), sample_sim_df[p[0]].dropna().tolist())
            print(p, s, pp)
        s = 'Mann-Whitney'
        add_stat_annotation(ax, data=sample_sim_df,
                    box_pairs=pairs, test=s)

        plt.show()

    ### num_clonotype / degree of sharing bar pplot
    if True:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace = 0.1)
        
        if target_exp == 'FMT1':
            if target_organ == 'LG':
                color = [color_list[i] for i in [0, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1]]
                df['degree of sharing'] = [90, 50, 40, 29, 18, 60, 70, 80, 100, 33, 22]
            elif target_organ == 'SP':
                df['degree of sharing'] = [100, 90, 80, 70, 60, 50, 40, 29, 18, 55, 44, 33, 22]
                color = [color_list[i] for i in [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,]]
        
        ax1.bar(x=df['degree of sharing'].tolist(), height=df['count'].tolist(), width=4, color=color, alpha=0.9)
        
        for p in ax1.patches: 
            if p.get_height() == 0:
                continue
            ax1.annotate("%d" % p.get_height(), (p.get_x() + p.get_width()/2., p.get_height() ), 
            ha='center', va='center', fontsize=16, color='black', xytext=(0, 10), 
            textcoords='offset points')
        ax1.set_ylabel("")
        ax1.set_xlabel("")
        
        ax2.bar(x=df['degree of sharing'].tolist(), height=df['count'].tolist(), width=4, color=color, alpha=0.9)
        for p in ax2.patches: 
            if p.get_height() == 0:
                continue
            ax2.annotate("%d" % p.get_height(), (p.get_x() + p.get_width()/2., p.get_height() ), 
            ha='center', va='center', fontsize=16, color='black', xytext=(0, 10), 
            textcoords='offset points')
        ax2.set_ylabel("")
        ax2.set_xlabel("Degree of sharing (%)")

        if target_exp == 'FMT1':
            if target_organ == 'LG':
                ylim1 = [120, 140, 160]
                ylim2 = [0, 20, 40]
            else:
                ylim1 = [200, 10000, 20000, 32000]
                ylim2 = [0, 50, 100]
        else:
            if target_organ == 'LG':
                ylim1 = [70, 80, 90]
                ylim2 = [0, 10, 20]
            else:
                ylim1 = [5000, 7000, 9000]
                ylim2 = [0, 300, 600]
        ax1.set_ylim(ylim1[0], ylim1[-1])
        ax2.set_ylim(ylim2[0], ylim2[-1])
        ax1.set_yticks(ylim1)
        ax2.set_yticks(ylim2)
        ax1.set_yticklabels(ylim1)
        ax2.set_yticklabels(ylim2)
        

        ax1.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax1.xaxis.tick_top()
        ax1.tick_params(labeltop=False)
        ax2.xaxis.tick_bottom()
        
        fig.text(0.01, 0.50, "Number of clonotypes", va='center', rotation = 'vertical')
        ax1.get_xaxis().set_visible(False)
        
        kwargs = dict(marker=[(-1, -0.7), (1, 0.7)], markersize=15,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
        ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
        
        plt.xticks([18, 22, 29, 33, 40, 50, 60, 70, 80, 90, 100], [20, 22, 30, 33, 40, 50, 60, 70, 80, 90, 100])

    ### frequency viloin plot   
    if True: 
        samples = ['NOD-FMT', 'PBS-T']
        
        fig, ax = plt.subplots(figsize=figsize)
        log_freq_dict = {k: [] for k in samples}
        area = {k: [] for k in samples}
        for s in samples:
            for i in range(len(freq_dict[s])):
                key = np.log10(freq_dict[s][i])
                degree = degree_dict[s][i]
                if degree < 2:
                    continue
                log_freq_dict[s].append(key)
                area[s].append(10*(5+key+degree*10)**2)
            
        plot_degree_dict = {s: [] for s in samples}
        for s in samples:
            for i in range(len(degree_dict[s])):
                degree = degree_dict[s][i]
                if degree < 2:
                    continue
                
                if s == samples[0]:
                    degree *= 10
                    if degree == 20:
                        degree = 18
                    elif degree == 30:
                        degree = 29
                    plot_degree_dict[s].append(degree-10)
                if s == samples[1]:
                    if degree == 2:
                        degree = 22
                    elif degree == 3:
                        degree = 33
                    elif degree == 4:
                        degree = 44
                    elif degree == 5:
                        degree = 55
                    else:
                        degree = degree * 100 / 9
                    plot_degree_dict[s].append(degree-10)
                    
        yticks = [-6, -5, -4, -3, -2, -1, 0]


        all_freq = np.array(log_freq_dict[samples[0]] + log_freq_dict[samples[1]])
        avg = np.mean(all_freq)
        std = np.std(all_freq)
        
        color = {}
        df_dict = {k: [] for k in ['sample', 'freq', 'degree']}
        for s in samples:
            for j in range(len(log_freq_dict[s])):
                df_dict['sample'].append(s)
                df_dict['freq'].append(log_freq_dict[s][j])
                degree = degree_dict[s][j]
                
                if s == samples[0]:
                    degree *= 10
                    if degree == 20:
                        degree = 18
                    elif degree == 30:
                        degree = 29
                    
                    df_dict['degree'].append(degree)
                else:
                    if degree == 2:
                        degree = 22
                    elif degree == 3:
                        degree = 33
                    elif degree == 4:
                        degree = 44
                    elif degree == 5:
                        degree = 55
                    else:
                        degree = degree *100 / 9
                    df_dict['degree'].append(degree)
        for i in range(10, 101):
            if i not in df_dict['degree']:
                df_dict['sample'].append(samples[-1])
                df_dict['freq'].append(-30)
                df_dict['degree'].append(i)
                
        df = pd.DataFrame(df_dict)

        color = [color_list[2]] * 8 + [color_list[0]] * 4 + [color_list[1]] * 7 + [color_list[0]] * 4 + [color_list[1]] * 7 + [color_list[0]] * 100
        color = [color_list[2]] * 8 + [color_list[0]] * 4 + [color_list[1]] * 7 + [color_list[0]] * 4 + [color_list[1]] * 7 + [color_list[0]] * 4 + [color_list[1]] * 6 +  [color_list[0]] * 4 + [color_list[1]] * 5 + [color_list[0]]*100
        v = sns.violinplot(x="degree", y="freq", data=df, ax=ax, legend=False, width=4, scale='width', palette=color)

        for violin in ax.collections:
            violin.set_alpha(0.2)
        ls = [l for i, l in enumerate(ax.lines)]
        [l.set_linestyle('solid') for i, l in enumerate(ls)] 
            
        ax2 = ax.twinx()

        ax2.axhline(avg, 0, 1, linestyle='--', color='gray')

        for i in range(len(samples)):
            s = samples[i]
            ax2.scatter(plot_degree_dict[s], log_freq_dict[s], c=color_list[i], alpha=0.7, s=50, label=s)
        
        fig.text(0.4, 0.4, "Avg frequency of all clonotypes", va='center', color='gray', fontsize=16)

        ax.set_xlim(2,96)
        ax2.set_xlim(2,96)
        ax.set_ylim(-6,0)
        ax2.set_ylim(-6,0)
        ax.set_ylabel("")
        ax2.set_ylabel("")
        ax2.set_xlabel("")
        ax2.yaxis.set_visible(False)
        
        ax.set_xlabel("Degree of sharing (%)")
  
        if target_organ == 'LG':
            plt.xticks([8, 12, 19, 23, 30, 40, 50, 60, 70, 80, 90], [20, 22, 30, 33, 40, 50, 60, 70, 80, 90, 100])
        else:
            plt.xticks([8, 12, 19, 23, 30, 34, 40, 45, 50, 60, 70, 80, 90], [20, 22, 30, 33, 40, 44, 50, 56, 60, 70, 80, 90, 100])
        plt.show()
        