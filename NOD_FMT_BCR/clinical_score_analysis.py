import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from statannot import add_stat_annotation
from itertools import combinations
from scipy.stats import ttest_ind, shapiro, mannwhitneyu, levene, friedmanchisquare, ttest_rel, normaltest, kstest
import numpy as np
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.formula.api as smf
import scikit_posthocs as sp

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 25

font = 'Arial'
figsize = (7, 7)

color_list = ['#B2182B', '#2166AC', '#808080']

def select_significance(list_1,list_2,list_3):
    if (shapiro(list_1)[1] > 0.05) and (shapiro(list_2)[1] > 0.05) and (shapiro(list_3)[1] > 0.05):
        if levene(list_1, list_2, list_3)[1] > 0.05:
            t,p = ttest_ind(list_1,list_2, equal_var=True)
            return 't-test_ind', p
        else:
            t,p = ttest_ind(list_1,list_2, equal_var=False)
            return 't-test_welch', p
    
    else:
        t,p = mannwhitneyu(list_1,list_2)
        return 'Mann-Whitney', p
        

if True:
    file = '/home/ajhml71531/sjogren_mouse/FMT/analysis/Opthalmology/NEI.csv'    
    df = pd.read_csv(file)
    groups = ['FMT', 'PBS']
    days = ['Pre-FMT', 'Post-FMT']
    final_plot_dict = {k: [] for k in ['Group', 'Value']}
    final_stat_dict = {k: [] for k in ['NOD-FMT', 'PBS-T']}
    for i in range(len(groups)):
        group = groups[i]
        stat_df = df[df['Sample'].str.contains(group)].iloc[:, 2:]
        stat_df.columns = ['Pre-FMT', '-', 'Post-FMT']
        stat_df = stat_df.drop(columns='-')
        plot_dict = {k: [] for k in ['Day', 'Value', 'Sample']}
        if i == 0:
            plot_dict['Day'].extend([days[0]]*20+[days[1]]*20)
            plot_dict['Sample'].extend(list(range(20))+list(range(20)))
            final_plot_dict['Group'].extend(['NOD-FMT']*20)
            final_stat_dict['NOD-FMT'].extend(list(stat_df[days[1]].dropna().values))
        elif i == 1:
            plot_dict['Day'].extend([days[0]]*18+[days[1]]*18)
            plot_dict['Sample'].extend(list(range(18))+list(range(18)))
            final_plot_dict['Group'].extend(['PBS-T']*18)
            final_stat_dict['PBS-T'].extend(list(stat_df[days[1]].dropna().values))
        final_plot_dict['Value'].extend(list(stat_df[days[1]].dropna().values))
        plot_dict['Value'].extend(list(stat_df[days[0]].values)+list(stat_df[days[1]].values))

        plot_df = pd.DataFrame(plot_dict)
        x_tick = []
        for d in plot_df['Day']:
            if d == days[0]:
                x_tick.append(1)
            elif d == days[1]:
                x_tick.append(2)
         
        plot_df['Day'] = x_tick
        
        fig, ax = plt.subplots(figsize=figsize)
    
        fig.text(-0.03, 0.50, "NEI score", va='center', rotation = 'vertical')
        
        mean = stat_df.mean()
        err = stat_df.sem()
        ax.errorbar(x=[0,1], y=mean.values, yerr=err.values, linestyle='None', capsize=20, color=color_list[i], elinewidth=3, capthick=3)

        sns.barplot(data=plot_df, x='Day', y='Value', color=color_list[i], alpha=0.8, ci=None)
        for patch in ax.patches :
            current_width = patch.get_width()
            diff = current_width - 0.5
            patch.set_width(0.5)
            patch.set_x(patch.get_x() + diff * .5)
        sns.stripplot(data=plot_df, x='Day', y='Value', size=7, color=color_list[i], ax=ax)
        plt.xlabel('')
        plt.ylabel('')

        ax.set_xticklabels(days)
        plt.yticks([0, 5,10,15])

if True:
    days = ['Pre-FMT', 'Post-FMT']
    PBS_before = [0.307402605,0.311045255,0.261979506,0.259787952,0.377347951, 0.356329865,0.323376773, 0.292930458,0.218333831, 0.286712451,0.235658043, 0.228408949,0.398092031, 0.340965208,0.300630734, 0.309461009,0.380447942, 0.31622276]#, 0.314826419, 0.362764757]
    PBS_after = [0.235395327,0.2308739,0.236016717, 0.2450504,0.278098185, 0.3074975,0.294614953, 0.2807517,0.245515165, 0.3241127,0.260554371, 0.3274627,0.355057887, 0.3568012,0.307121243, 0.3132943,0.317132008, 0.3271834]
    
    FMT_before = [0.449113924, 0.271898734,0.390638298, 0.30443769,0.45245509, 0.346886228,0.32656464, 0.390817603,0.396050574, 0.360466632,0.345714286, 0.364464286,0.403229974, 0.356459948,0.301903629, 0.38958953,0.328119763, 0.284102335,0.480955297, 0.489895897]
    FMT_after = [0.3122785, 0.2980079,0.220068, 0.2782328,0.3026835, 0.2105148,0.363159, 0.31408,0.3196255, 0.3035157,0.3810017, 0.400261,0.5003829, 0.4400817,0.2844012, 0.3600746,0.3259941, 0.2464248,0.2582305, 0.2976992]
    
    NC_before = [0.401161665,0.311907067,0.282591631, 0.333313206,0.296103896, 0.283623693,0.357701047, 0.312863714,0.343426791, 0.383364486]
    NC_after = [0.293771817,0.3428716,0.350633549, 0.3491231,0.282552579, 0.3010029,0.34090723, 0.3182265,0.404247692, 0.350345]

    plot_dict = {k: [] for k in ['Group', 'Day', 'Value']}
    plot_dict['Group'].extend(['FMT']*40+['PBS']*36)
    plot_dict['Day'].extend([days[0]]*20+[days[1]]*20+[days[0]]*18+[days[1]]*18)
    plot_dict['Value'].extend(FMT_before+FMT_after+PBS_before+PBS_after)

    plot_df = pd.DataFrame(plot_dict)
    
    groups = ['FMT', 'PBS']
    for i in range(len(groups)):
        group = groups[i]
        p_df = plot_df[plot_df['Group'] == group]
        
        fig, ax = plt.subplots(figsize=figsize)
    
        fig.text(-0.03, 0.50, "Tear production/body weight", va='center', rotation = 'vertical')
   
        s = 't-test_paired'
        if i == 0:
            stat_df = pd.DataFrame({days[0]: FMT_before, days[1]: FMT_after})
        elif i == 1:
            stat_df = pd.DataFrame({days[0]: PBS_before, days[1]: PBS_after})
        
        p_df['Day'] = [2 if d == days[0] else 3 for d in p_df['Day'].tolist()]

        mean = stat_df.mean()
        err = stat_df.sem()
        ax.errorbar(x=[0,1], y=mean.values, yerr=err.values, linestyle='None', capsize=20, color=color_list[i], elinewidth=3, capthick=3)
        sns.barplot(data=p_df, x='Day', y='Value', color=color_list[i], alpha=0.8, ci=None)
        for patch in ax.patches :
            current_width = patch.get_width()
            diff = current_width - 0.5
            patch.set_width(0.5)
            patch.set_x(patch.get_x() + diff * .5)
        sns.stripplot(data=p_df, x='Day', y='Value', size=7, color=color_list[i], ax=ax)
        plt.xlabel('')
        plt.ylabel('')

        ax.set_xticklabels(days)
        plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])


if True:
    IL6_PBS = [0.099580272,0.528595998,0.482022748,0.985692099,1.597019226,1.284948527,1.953306792,1.197082397,0.871751941]
    IL6_FMT = [3.005365897,2.511324252,0.335038895,1.692661601,2.92988729,9.365567834,5.274999276,3.414716317,15.5906515,1.244297515]

    Muc_PBS = [0.657917385,0.330168214,0.888805845,0.816327233,0.965276444,0.963945045,2.052691701,1.435459273,0.889408859]
    Muc_FMT = [0.130220287,0.039546515,0.119438108,0.23044271,0.081258765,0.711640955,0.607286934,0.1642679,0.586380613,0.006123654]

    plot_dict = {k: [] for k in ['Group', 'Type', 'Value']}
    plot_dict['Group'].extend(['PBS']*9+['FMT']*10+['PBS']*9+['FMT']*10)
    plot_dict['Type'].extend([0]*19+[1]*19)
    plot_dict['Value'].extend(IL6_PBS+IL6_FMT+Muc_PBS+Muc_FMT)

    plot_df = pd.DataFrame(plot_dict)
   
    for i in range(2):
        fig, ax = plt.subplots(figsize=figsize)

        fig.text(-0.03, 0.50, "mRNA level (Relative units)", va='center', rotation = 'vertical')
        s = 'Mann-Whitney'
        if i == 0:
            stat_df = pd.DataFrame.from_dict({'NOD-FMT': IL6_FMT, 'PBS-T': IL6_PBS}, orient='index').T
        else:
            stat_df = pd.DataFrame.from_dict({'NOD-FMT': Muc_FMT, 'PBS-T': Muc_PBS}, orient='index').T
        
        mean = stat_df.mean()
        err = stat_df.sem()
        ax.errorbar(x=[0], y=mean.values[0], yerr=err.values[0], linestyle='None', capsize=20, color=color_list[0], elinewidth=3, capthick=3)
        ax.errorbar(x=[1], y=mean.values[1], yerr=err.values[1], linestyle='None', capsize=20, color=color_list[1], elinewidth=3, capthick=3)
       
        p_df = plot_df[plot_df['Type'] == i]
        p_df['Group'] = [2 if d == 'FMT' else 3 for d in p_df['Group'].tolist()]
        sns.barplot(data=p_df, x='Group', y='Value', palette={2: color_list[0], 3: color_list[1]}, alpha=0.8, ci=None)
        for patch in ax.patches :
            current_width = patch.get_width()
            diff = current_width - 0.5
            patch.set_width(0.5)
            patch.set_x(patch.get_x() + diff * .5)
        sns.stripplot(data=p_df, x='Group', y='Value', size=7, palette={2: color_list[0], 3: color_list[1]}, ax=ax)
        plt.xlabel('')
        plt.ylabel('')

        # ax.set_xticks([1,3])
        ax.set_xticklabels(['NOD-FMT', 'PBS-T'])
        add_stat_annotation(ax, data=stat_df,
                        box_pairs=[('NOD-FMT', 'PBS-T')], test=s)
        if i == 0:
            plt.yticks([0, 5,10,15,20])
        else:
            plt.yticks([0, 0.5,1,1.5,2,2.5])