# coding=utf-8
import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State


logo = "https://raw.githubusercontent.com/sradiouy/UTRme/master/Images/Logo2.png"




# row1
header = html.Div([
    html.Img(className="logo",src=logo),
], className='navbar navbar-default navbar-static-top',style={"background-color": "white"})



configuration = html.Div([
                html.Div([
                    html.Br([]),
                    html.H6('UTRme: A Scoring-Based Tool to Annotate Untranslated Regions in Trypanosomatid Genomes',
                            className="gs-header gs-table-header padded",style={"margin-top":"0px","text-align": "center","font-size":"17px"}),
                    html.Br([]),
                    html.P("\
                           Most signals involved in post-transcriptional regulatory networks are located in the untranslated regions (UTRs) of the mRNAs. Therefore, to deepen our understanding of gene expression regulation, delimitation of these regions with high accuracy is needed. In this context, the definition of the UTR regions becomes of key importance.UTRme implements a multiple scoring system tailored to address the issue of false positive UTR assignment that frequently arise because of the characteristics of the intergenic regions. Even though it was developed for trypanosomatids, the tool can be used to predict 3′ sites in any eukaryote and 5′ UTRs in any organism where trans-splicing occurs (such as the model organism C. elegans)."),
                ], className="twelve columns about_text"),
                html.Div([
                    html.H6([""],
                            className="gs-header gs-table-header padded",style={"height": "22px"}),
                ], className="twelve columns about_text"),
            ], className="row ")

layout_home = html.Div([
        header,
        configuration])