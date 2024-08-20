#!/usr/bin/env python
"""
Plotly dash application to visualize some results from step2

usage:
python app_step2.py -i [/path/to/directory/cluster_analysis, required] -k [keggBrite, optional] -c [cogFunction, optional]
"""
import dash
from dash import Dash, html, dcc, callback, Output, Input, Patch, ctx
from dash_bio import Clustergram
import dash_bootstrap_components as dbc
import dash_ag_grid as dag
import polars as pl
import plotly.express as px
import plotly.graph_objects as go
import io
from math import ceil
import argparse
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--inputDir", required=True)
ap.add_argument("-k", "--keggBrite", required=False, default="KeggBrite.tsv")
ap.add_argument("-c", "--cogFunction", required=False, default="cogFunction.tsv")
args = vars(ap.parse_args())

# Data importation
folder_path = Path(args["inputDir"])
if not folder_path.is_dir():
    raise NotADirectoryError("The given argument --inputDir (-i) must be a directory")

clusters_info = pl.scan_csv(f"{folder_path}/clusters_info.tsv", separator='\t')
clusters_abundance = pl.scan_csv(f"{folder_path}/clusters_abundance.tsv", separator='\t', infer_schema_length=0)
clusters_annotation = pl.scan_csv(f"{folder_path}/clusters_annotation.tsv", separator='\t', infer_schema_length=0)
bin_annotation = pl.scan_csv(f"{folder_path}/bin_annotation_all.tsv", separator='\t')
contigs2bins = pl.scan_csv(f"{folder_path}/contigs2bins_all.tsv", separator='\t')
keggBritePath = Path(args['keggBrite'])
if keggBritePath.is_file():
    keggBrite = pl.scan_csv(keggBritePath, separator='\t')
cogFuncPath = Path(args["cogFunction"])
if cogFuncPath.is_file():
    cogFunction = pl.scan_csv(cogFuncPath, separator='\t')

clusters_abundance = clusters_abundance.with_columns(pl.col("*").exclude("cluster").cast(pl.Int32))
clusters_annotation = clusters_annotation.with_columns(pl.col("^.*_E_value$").cast(pl.Float64))
annot_names = pl.from_dict({"colnames":clusters_annotation.collect_schema().names()}).to_series()
if annot_names.str.contains("Kegg_id").any() and keggBritePath.is_file():
    clusters_annotation = clusters_annotation.join(keggBrite, left_on="Kegg_id", right_on="keggId", coalesce=True, how="left")
    clusters_annotation = clusters_annotation.rename({"briteName":"Kegg_Brite_name"})
if annot_names.str.contains("cog-20_name").any() and cogFuncPath.is_file():
    clusters_annotation = clusters_annotation.join(cogFunction, left_on="cog-20_name", right_on="cog_name", coalesce=True, how="left")

samples = clusters_info.select(pl.col('source')).unique().collect().to_series().sort().to_list()

annot_names = pl.from_dict({"colnames":clusters_annotation.collect_schema().names()}).to_series()
annot_names = annot_names.filter(annot_names.str.ends_with("_name")).str.strip_suffix("_name").sort().to_list()

sum_expr = sum(pl.col(col) for col in samples)

taxon_1_list = [{"label":"Domain", "value":"d__"},
                 {"label":"Phylum", "value":"p__"},
                 {"label":"Class", "value":"c__"},
                 {"label":"Order", "value":"o__"},
                 {"label":"Family", "value":"f__"},
                 {"label":"Genus", "value":"g__"},
                 {"label":"Species", "value":"s__"}]

fig = go.Figure(layout=dict(template='plotly'))


df_abundance = clusters_abundance.fill_null(0).with_columns([sum_expr.alias("count")]).sort(by="count", descending=True).select(pl.col(samples), pl.col("cluster")).join(clusters_annotation.select(["cluster", "cog_func_name"]), on="cluster", how="left").drop("cluster").group_by(by=pl.col("cog_func_name")).sum().drop("cog_func_name").rename({"by":"cog_func_name"}).collect().to_pandas().set_index("cog_func_name")
columns = list(df_abundance.columns.values)
rows = list(df_abundance.index)
fig_abund = Clustergram(data=df_abundance.loc[rows].values, row_labels=rows, column_labels=columns, height=800, width=700, line_width=0.7, center_values=False)

# test_df = clusters_info.select(["seq_id", "cluster", "contigId", pl.col("seq_name_fasta").str.extract(r";gc_cont=(.*$)", group_index=1).cast(pl.Float64).alias("gc_cont")]).join(clusters_abundance.fill_null(0).with_columns([sum_expr.alias("count")]), on="cluster", coalesce=True).join(contigs2bins, on="contigId", coalesce=True).join(bin_annotation.with_columns("binId", pl.col('classification').str.extract(r"o__\w+|^Unclassified \w+$", group_index=0).str.strip_prefix("o__").fill_null("None")), on="binId", coalesce=True).select(["cluster", "gc_cont", "count", "classification"]).filter(pl.col("count") != 0).collect()
# fig_abund = px.scatter(test_df, x="gc_cont", y="count", color="classification", opacity=0.5)


app = Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.CERULEAN], pages_folder="")


# Graphs page layout:
dash.register_page("Graphs", path='/', layout=html.Div([
    html.Div([
        dcc.Markdown("**Clusters Info**"),
    ]),
    html.Hr(),
    html.Div([
        "Bin characteristic:",
        dcc.Dropdown(['Completeness', 'Contamination', 'num_seqs', 'N50', 'sum_len'], 'Completeness', id='yaxis'),
        "Sample:",
        dcc.Dropdown(samples + ['all'], 'all', id='sample'),
        dcc.RadioItems(['GTDB-Tk classification', 'Checkm lineage'], 'GTDB-Tk classification', id='classif_type', inline=False),
        dcc.Dropdown(id='classif_list'),
        dcc.Graph(id='graph_bin'),
        "Boxplot x-axis:",
        dcc.RadioItems(['Sample', 'Classification'], 'Sample', id='xaxis_box', inline=False),
        dcc.Store(id='classif_lf', storage_type='memory'),
        dcc.Graph(id='box_bin')
    ]),
    html.Hr(),
    html.Div([
        "Type of annotation:",
        dcc.Dropdown(annot_names, annot_names[0], id='select_annot'),
        "Sample:",
        dcc.Dropdown(samples + ['all'], 'all', id='sample_annot'),
        dcc.RadioItems(['Abundance', 'Number of genes'], 'Abundance', id='count_type', inline=False),
        dcc.Graph(id='graph_annot'),
        dcc.Markdown(id='slider_name'),
        dcc.RangeSlider(step = 1, id='count_slider',
                         tooltip={"placement": "bottom", "always_visible": True}),
        dcc.Store(id="lf_store", storage_type='memory'),
    ]),
    html.Hr(),
    html.Div([
        dcc.Graph(figure=fig_abund, id='graph_test')
    ]),
], style={'padding': 10}))


# Download page layout
dash.register_page("download", layout=html.Div([
    dcc.Markdown("**Download genes**"),
    html.Div([
        html.Div([
            "Taxon:",
            html.Div([
                html.Div([
                    dcc.Dropdown(taxon_1_list, "g__", searchable=False, id="taxon_1"),
                ], style={'padding': 0, 'flex': 1}),
                html.Div([
                    dcc.Dropdown(id="taxon_2"),
                    dcc.Store(id='taxon_2_lf', storage_type='memory'),
                ], style={'padding': 0, 'flex': 1}),
            ], style={'display': 'flex', 'flexDirection': 'row'}),
            "Sample:",
            dcc.Dropdown(value='all', id="sample"),
            "Bin:",
            dcc.Dropdown(value='all', id="bin"),
        ], style={'padding': 10, 'flex': 1}),
        html.Div([
            dcc.Checklist(options=[{"label":" Only annotated genes", "value":"only_annot"}], value=["only_annot"], id="checklist_genes"),
            html.Div(id='test-output-slider'),
            dcc.Slider(1, 54, 0.001, value=1, marks={i: '{}'.format(round(1/(i * 10**i), round(i)+2)) for i in range(1, 55, 5)}, id="E_value_slider", disabled=False),
        ], style={'padding': 10, 'flex': 1}),
    ], style={'display': 'flex', 'flexDirection': 'row'}),
    html.Br(),
    dcc.Input(id="input-row", placeholder="Quick filter..."),
    dcc.Loading(
        html.Div([
            dag.AgGrid(
                id="row-selection",
                columnSize="sizeToFit",
                defaultColDef={"filter": True},
                #getRowId="params.data.id",
                dashGridOptions={
                    "rowSelection": "multiple",
                    "pagination": True,
                    "paginationAutoPageSize": True,
                    "animateRows": False,
                },
            ),
        ]),
        target_components ={"row-selection": "*"}
    ),
    html.Br(),
    html.Button("Download FASTA", id="btn_fasta"),
    dcc.Download(id="download_fasta"),
    dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("No genes selected")),
                dbc.ModalBody("Please select genes to download"),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close", id="close", className="ms-auto", n_clicks=0
                    )
                ),
            ],
            id="no_download_window",
            is_open=False,
        ),
], style={'padding': 10}))


# Home page layout
app.layout = html.Div([
    html.Div([
        dbc.NavbarSimple(
        children=[
            dbc.NavItem(dbc.NavLink(f"{page['name']}", href=page["relative_path"])) for page in dash.page_registry.values()
        ],
        brand="Clustermg Step 2 Results",
        brand_href="#",
        color="primary",
        dark=True,
        ),
    ]),
    dash.page_container
])


# Graphs page callbacks:
@callback(
        Output(component_id='classif_list', component_property='options'),
        Output(component_id='classif_list', component_property='value'),
        Output(component_id='classif_lf', component_property='data'),
        Input(component_id='classif_type', component_property='value')
)
def graph_update_classif_type(classif_type):
    if classif_type == 'GTDB-Tk classification':
        bin_annot_lf = bin_annotation.with_columns(pl.col('classification').str.extract(r"g__\w+|^Unclassified \w+$", group_index=0).str.strip_prefix("g__").fill_null("None"))
        options = bin_annot_lf.select(pl.col('classification').unique()).collect().to_series().to_list()
    elif classif_type == 'Checkm lineage':
        bin_annot_lf = bin_annotation.with_columns(pl.col('checkm_lineage').str.extract(r"[a-zA-Z]__\w+", group_index=0).fill_null("None"))
        options = bin_annot_lf.select(pl.col('checkm_lineage').unique()).collect().to_series().to_list()
    options = options + ['all']
    value = 'all'
    return options, value, bin_annot_lf.serialize(format='json')

@callback(
    Output(component_id='graph_bin', component_property='figure'),
    Input(component_id='yaxis', component_property='value'),
    Input(component_id='sample', component_property='value'),
    Input(component_id='classif_lf', component_property='data'),
    Input(component_id='classif_type', component_property='value'),
    Input(component_id='classif_list', component_property='value')
)
def graph_update_bar_bin(col_chosen, sample_name, classif_lf, classif_type, classif_list):
    if sample_name != 'all' and sample_name != None:
        bin_annotation_filt = pl.LazyFrame.deserialize(io.StringIO(classif_lf), format='json').filter(pl.col('binId').str.starts_with(sample_name + '|'))
    else:
        bin_annotation_filt = pl.LazyFrame.deserialize(io.StringIO(classif_lf), format='json')
    if classif_type == 'GTDB-Tk classification':
        color_name = 'classification'
        legend_name = 'Genus'
    elif classif_type == 'Checkm lineage':
        color_name = 'checkm_lineage'
        legend_name = 'Lineage'
    if classif_list != 'all' and classif_list != None:
        bin_annotation_filt = bin_annotation_filt.filter(pl.col(color_name).str.contains(classif_list, literal=True))
    bin_annot_df = bin_annotation_filt.collect()
    list_bin = bin_annot_df['binId'].unique().sort().to_list()
    fig = px.bar(bin_annot_df, x='binId', y=col_chosen, color=color_name, category_orders={"binId":list_bin})
    fig.update_layout(legend_title_text=legend_name)
    return fig

@callback(
    Output(component_id='box_bin', component_property='figure'),
    Input(component_id='yaxis', component_property='value'),
    Input(component_id='xaxis_box', component_property='value'),
    Input(component_id='sample', component_property='value'),
    Input(component_id='classif_list', component_property='value'),
    Input(component_id='classif_type', component_property='value'),
)
def graph_update_box_bin(col_chosen, xaxis, sample_name, classif, classif_type):
    if classif_type == 'GTDB-Tk classification':
        bin_annot_class = bin_annotation.with_columns(pl.col('classification').str.extract(r"g__\w+|^Unclassified \w+$", group_index=0).str.strip_prefix("g__").fill_null("None"))
    elif classif_type == 'Checkm lineage':
        bin_annot_class = bin_annotation.with_columns(pl.col('checkm_lineage').str.extract(r"[a-zA-Z]__\w+", group_index=0).fill_null("None").alias('classification'))
    if xaxis == 'Sample':
        if sample_name != 'all' and sample_name != None:
            bin_annot_class = bin_annot_class.filter(pl.col('binId').str.starts_with(sample_name + '|'))
        df_bin = bin_annot_class.join(contigs2bins, on='binId', how='left', coalesce=True).join(clusters_info.select(pl.col('contigId', 'source')), on='contigId', how='left', coalesce=True).unique(subset=["binId", "source"], maintain_order=True).collect()
        x_name = 'source'
        x_list = samples
        color_name = 'classification'
    elif xaxis == 'Classification':
        if classif != 'all' and classif != None:
            bin_annot_class = bin_annot_class.filter(pl.col('classification').str.starts_with(classif))
        df_bin = bin_annot_class.join(contigs2bins, on='binId', how='left', coalesce=True).join(clusters_info.select(pl.col('contigId', 'source')), on='contigId', how='left', coalesce=True).unique(subset=["binId", "source"], maintain_order=True).collect()
        x_name = 'classification'
        x_list = df_bin['classification'].unique().to_list()
        color_name = 'source'

    list_bin = df_bin[x_name].unique().sort().to_list()

    fig = px.strip(df_bin,
         x=x_name,
         y=col_chosen,
         color=color_name,
         stripmode='overlay',
         category_orders={x_name:list_bin})
    fig.update_traces({'marker':{'size': 8}})
    for x in x_list:
        fig.add_trace(go.Box(y=df_bin.filter(pl.col(x_name) == x)[col_chosen], name=x, opacity=0.6, boxpoints=False, jitter=0.4, pointpos=-2))
    fig.update_layout(legend={'traceorder':'normal'})
    return fig

@callback(
        Output(component_id='count_slider', component_property='min'),
        Output(component_id='count_slider', component_property='max'),
        Output(component_id='count_slider', component_property='marks'),
        Output(component_id='count_slider', component_property='value'),
        Output(component_id='lf_store', component_property='data'),
        Output(component_id='slider_name', component_property='children'),
        Input(component_id='select_annot', component_property='value'),
        Input(component_id='count_type', component_property='value'),
        Input(component_id='sample_annot', component_property='value')
)
def graph_update_annot_type(select_annot, count_type, sample_annot):
    if count_type == "Number of genes":
        if sample_annot != 'all' and sample_annot != None:
                annot_count = clusters_abundance.select(pl.col("cluster"), pl.col(sample_annot).is_not_null()).join(clusters_annotation, on="cluster").filter(pl.col(sample_annot)).select(pl.col(f"{select_annot}_name")).drop_nulls().group_by(f"{select_annot}_name").len('count')
        else:
            test_lf = clusters_abundance.select(pl.col("cluster"), pl.col(samples).is_not_null()).join(clusters_annotation, on="cluster")
            final_lf = pl.LazyFrame.clear(test_lf.select(pl.col(f"{select_annot}_name")))
            for col in samples:
                new_lf = test_lf.rename({col: "sample"}).filter(pl.col("sample")).group_by([f"{select_annot}_name", "sample"]).len(col).drop("sample").cast({pl.UInt32: pl.Int32})
                final_lf = final_lf.join(new_lf, on=f"{select_annot}_name", how="full", coalesce=True)
            annot_count = final_lf.fill_null(0).with_columns([sum_expr.alias("count")]).drop_nulls()
    elif count_type == "Abundance":
        if sample_annot != 'all' and sample_annot != None:
            annot_count = clusters_abundance.select(pl.col(["cluster", sample_annot])).drop_nulls().join(clusters_annotation, on="cluster").select([sample_annot, f"{select_annot}_name"]).drop_nulls().group_by(f"{select_annot}_name").agg(pl.col(sample_annot).sum().alias("count"))
        else:
            annot_count = clusters_abundance.join(clusters_annotation, on="cluster").group_by(f"{select_annot}_name").agg(pl.col(samples).sum()).with_columns([sum_expr.alias("count")]).drop_nulls()
    annot_count = annot_count.filter(pl.col("count") > 0)
    annot_count = annot_count.sort(by="count", descending=True)
    min = annot_count.select(pl.col('count').min()).collect()[0, 0]
    max = annot_count.select(pl.col('count').max()).collect()[0, 0]
    marks = {i: str(i) for i in range(0, max, 10 * ceil((max/15)/10))}
    value = [min, max]
    children = f"Min/Max {count_type}:"
    return min, max, marks, value, annot_count.serialize(format='json'), children

@callback(
    Output(component_id='graph_annot', component_property='figure'),
    Input(component_id='count_slider', component_property='value'),
    Input(component_id='select_annot', component_property='value'),
    Input(component_id='lf_store', component_property='data'),
    Input(component_id='count_type', component_property='value'),
    Input(component_id='sample_annot', component_property='value')
)
def garph_update_annot_hist(count_range, selected_type, count_lf, count_type, sample_annot):
    annot_df = pl.LazyFrame.deserialize(io.StringIO(count_lf), format='json').collect()
    annot_df = annot_df.filter((pl.col("count") >= count_range[0]) & (pl.col("count") <= count_range[1]))
    if sample_annot != 'all' and sample_annot != None:
        fig = px.histogram(annot_df, x=f"{selected_type}_name", y="count")
    else:
        fig = px.histogram(annot_df, x=f"{selected_type}_name", y=samples)
    fig.update_layout(yaxis_title=count_type, legend_title_text="Sample")
    return fig


# Download page callbacks
@callback(
        Output("taxon_2", "options"),
        Output("taxon_2_lf", "data"),
        Input("taxon_1", "value")
)
def downl_update_taxon_2(taxon_1):
    bin_annot_class = bin_annotation.with_columns(pl.col('classification').str.extract(fr"{taxon_1}\w+|^Unclassified \w+$", group_index=0).str.strip_prefix(f"{taxon_1}").fill_null("None"))
    return bin_annot_class.select(pl.col("classification").unique()).collect().to_series().to_list(), bin_annot_class.serialize(format='json')

@callback(
        Output("sample", "options"),
        Input("taxon_2", "value"),
        Input("taxon_2_lf", "data")
)
def downl_update_sample(taxon_2, taxon_2_lf):
    if taxon_2 == None:
        options = samples + ['all']
    else:
        bin_annotation_filt = pl.LazyFrame.deserialize(io.StringIO(taxon_2_lf), format='json').filter(pl.col('classification').str.starts_with(taxon_2))
        bin_annotation_filt.select(pl.col("binId").str.extract(r"(.*)\|", group_index=1)).unique().collect().to_series().to_list()
        options = bin_annotation_filt.select(pl.col("binId").str.extract(r"(.*)\|", group_index=1)).unique().collect().to_series().to_list() + ['all']
    return options

@callback(
        Output("bin", "options"),
        Input("taxon_2", "value"),
        Input("sample", "value"),
        Input("taxon_2_lf", "data")
)
def downl_update_bin(taxon_2, sample, taxon_2_lf):
    bin_annotation_filt = bin_annotation
    if taxon_2 != None:
        bin_annotation_filt = pl.LazyFrame.deserialize(io.StringIO(taxon_2_lf), format='json').filter(pl.col('classification').str.starts_with(taxon_2))
    if sample != 'all' and sample != None:
        bin_annotation_filt = bin_annotation_filt.filter(pl.col("binId").str.contains(fr"({sample})\|"))
    return bin_annotation_filt.select(pl.col("binId")).collect().to_series().to_list() + ["all"]

@callback(
        Output("E_value_slider", "disabled"),
        Input("checklist_genes", "value")
)
def downl_disable_slider(check):
    if check:
        return False
    else:
        return True

@callback(
        Output('test-output-slider', 'children'),
        Input('E_value_slider', 'value')
)
def downl_display_e_value(value):
    return f"Min E-value: {round(1/(value * 10**value), round(value)+2)}"

@callback(
    Output("row-selection", "dashGridOptions"),
    Output("row-selection", "columnDefs"),
    Output("row-selection", "rowData"),
    Input("taxon_2", "value"),
    Input("sample", "value"),
    Input("bin", "value"),
    Input("E_value_slider", "value"),
    Input("E_value_slider", "disabled"),
    Input("input-row", "value"),
    Input("taxon_2_lf", "data")
)
def downl_update_filter(taxon_2, sample, bin, min_e_value, slider_is_dis, filter_value, taxon_2_lf):
    selection_lf = pl.LazyFrame.deserialize(io.StringIO(taxon_2_lf), format='json').select(["binId", "classification"]).join(contigs2bins, on='binId', coalesce=True).join(clusters_info.select(["seq_id", "contigId"]), on='contigId', coalesce=True).join(clusters_annotation, on='seq_id', coalesce=True)
    if not slider_is_dis:
        min_e_value = round(1/(min_e_value * 10**min_e_value), round(min_e_value)+2)
        sele_e_value = selection_lf.select(pl.col(pl.Float64) < min_e_value).with_columns(any=pl.any_horizontal("*")).select(pl.col("any")).collect().to_series()
        selection_lf = selection_lf.filter(sele_e_value)
    if taxon_2 != None:
        selection_lf = selection_lf.filter(pl.col("classification").str.contains(taxon_2, literal=True))
    if sample != 'all' and sample != None:
        selection_lf = selection_lf.filter(pl.col("binId").str.contains(fr"({sample})\|"))
    if bin != 'all' and bin != None:
        selection_lf = selection_lf.filter(pl.col("binId").str.contains(bin, literal=True))
    columnDefs = [{"field": "seq_id", "checkboxSelection": True, "headerCheckboxSelection": True}]
    for col in selection_lf.collect_schema().names():
        if ("_name" in col or "_id" in col or "_E_value" in col) and "seq_id" not in col:
            columnDefs += [{"field": col, "checkboxSelection": False}]
    gridOptions_patch = Patch()
    gridOptions_patch["quickFilterText"] = filter_value
    return gridOptions_patch, columnDefs, selection_lf.select(pl.selectors.matches(r"_(id|name|E_value)$")).collect().to_dicts()

@callback(
    Output("download_fasta", "data"),
    Output("no_download_window", "is_open"),
    Input("btn_fasta", "n_clicks"),
    Input("row-selection", "selectedRows"),
    prevent_initial_call=True,
)
def downl_download_fasta(n_clicks, selectedRows):
    if "btn_fasta" == ctx.triggered_id:
        if selectedRows != []:
            selected_ids = [row['seq_id'] for row in selectedRows]
            selected_seqs = pl.LazyFrame({"seq_id":selected_ids}).join(clusters_info, on="seq_id", coalesce=True).select(["seq_name_fasta", "seq"])
            fasta = ""
            for row in selected_seqs.collect().rows():
                fasta = fasta + f"{row[0]}\n{row[1]}\n"
            return dict(content=fasta, filename="selected_genes.fasta"), False
        else:
            return None, True
    else:
        return None, False

@callback(
    Output("no_download_window", "is_open", allow_duplicate=True),
    Input("close", "n_clicks"),
    prevent_initial_call=True
)
def downl_close_window(value):
    return False


if __name__ == '__main__':
    app.run(debug=True)