import numpy as np
import pandas as pd
import seaborn as sns


def compute_qc(adata):
    import scanpy as sc
    print('Calculate QC stats...')
    if 'feature_name' in adata.var.columns:
        var_names = adata.var['feature_name']
    else:
        var_names = adata.var_names

    adata.var["mito"] = var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)


def plot_qc_joint(
    df,
    x,
    y,
    log_x=1,
    log_y=1,
    hue=None,
    main_plot_function=None,
    marginal_hue=None,
    marginal_legend=False,
    x_threshold=None,
    y_threshold=None,
    title='',
    return_df=False,
    **kwargs,
):
    """
    Plot scatter plot with marginal histograms from df columns.

    :param df: observation dataframe
    :param x: df column for x axis
    :param y: df column for y axis
    :param log: log base for transforming values. Default 1, no transformation
    :param hue: df column with annotations for color coding scatter plot points
    :param marginal_hue: df column with annotations for color coding marginal plot distributions
    :param x_threshold: tuple of upper and lower filter thresholds for x axis
    :param y_threshold: tuple of upper and lower filter thresholds for y axis
    :param title: Title text for plot
    :return:
        seaborn plot (and df dataframe with updated values, if `return_df=True`)
    """
    if main_plot_function is None:
        main_plot_function = sns.scatterplot
    if pd.api.types.is_numeric_dtype(df[hue]):
        kwargs |= {'palette': 'plasma', 'legend': 'brief'}
    else:
        kwargs |= {'palette': None, 'legend': df[hue].nunique() <= 20}
    if not x_threshold:
        x_threshold=(0, np.inf)
    if not y_threshold:
        y_threshold=(0, np.inf)

    def log1p_base(_x, base):
        return np.log1p(_x) / np.log(base)

    if log_x > 1:
        x_log = f'log{log_x} {x}'
        df[x_log] = log1p_base(df[x], log_x)
        x_threshold = log1p_base(x_threshold, log_x)
        x = x_log
    
    if log_y > 1:
        y_log = f'log{log_y} {y}'
        df[y_log] = log1p_base(df[y], log_y)
        y_threshold = log1p_base(y_threshold, log_y)
        y = y_log

    g = sns.JointGrid(
        data=df,
        x=x,
        y=y,
        xlim=(0, df[x].max()),
        ylim=(0, df[y].max()),
    )
    # main plot
    g.plot_joint(
        main_plot_function,
        data=df,
        hue=hue,
        **kwargs,
    )
    
    # marginal hist plot
    if marginal_hue in df.columns:
        marginal_hue = None if df[marginal_hue].nunique() > 100 else marginal_hue
    use_marg_hue = marginal_hue is not None
    g.plot_marginals(
        sns.histplot,
        data=df,
        hue=marginal_hue,
        legend=marginal_legend,
        element='step' if use_marg_hue else 'bars',
        fill=False,
        bins=100
    )

    g.fig.suptitle(title, fontsize=12)

    # x threshold
    for t, t_def in zip(x_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axvline(x=t, color='red')
            g.ax_marg_x.axvline(x=t, color='red')

    # y threshold
    for t, t_def in zip(y_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axhline(y=t, color='red')
            g.ax_marg_y.axhline(y=t, color='red')

    if return_df:
        return g, df
    return g
