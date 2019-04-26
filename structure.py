import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import os

import header
from app import app,server


def check_adapter_value(adapter):
    if adapter != "":
        valid_dna = 'ACGT'
        sequence = adapter.upper()
        return all(i in valid_dna for i in sequence)
    else:
        return True

def check_core(core):
    if core == "":
        return False
    try:
        val = int(core)
    except ValueError:
        return False
    return True

def check_id(id_value):
    if id_value == "":
        return False
    return True

def check_basename(basename_value):
    if basename_value[0].isdigit():
        return basename_value,False
    if basename_value == "":
        basename_value = "UTRme-Run"
    return basename_value,True

def check_sl_value(sl_value,organism_value):
    sl_line = ""
    if organism_value != "Other" and sl_value != "":
        sl_line = "You have selected ",organism_value," and specify a SL sequence. ", "If you want to input a different SL, please select Other."
        return sl_line,False
    elif organism_value == "Other" and not check_adapter_value(sl_value):
        sl_line = "Your SL sequence is not valid"
        return sl_line,False
    elif organism_value == "Other" and sl_value == "":
        sl_line = "Your must write a SL sequence if you select Other"
        return sl_line,False
    else:
        return sl_line,True

def check_annotation(annotation_value):
    annotation_line = ""
    if annotation_value == "GFF format":
        annotation_line = "You must put a valid path. Change GFF format"
        return annotation_line, False
    elif annotation_value == "":
        annotation_line = "You must put a annotation valid path."
        return annotation_line, False
    else:
        return annotation_line, True

def check_genome(genome_value):
    genome_line = ""
    if genome_value == "Fasta format":
        genome_line = "You must put a valid path. Change Fasta format"
        return genome_line, False
    elif genome_value == "":
        genome_line = "You must put a genome valid path."
        return genome_line, False
    else:
        return genome_line, True

def check_second(second_value,experiment_value):
    second_line = ""
    if experiment_value == "Paired-end":
        if second_value == "Full path to the folder where the fastq files (gzipped or not) are located. Second Pair or leave empty if single-end" or second_value == "":
            second_line = "You must put a valid path to the folder where the fastq files of the second-pair are located."
            return second_line,second_value,False
        else:
            return second_line,second_value,True
    else:
        if second_value == "Full path to the folder where the fastq files (gzipped or not) are located. Second Pair or leave empty if single-end" or second_value == "":
            second_value = ""
            return second_line,second_value,True
        else:
            second_line = "You must leave empty Second-pair folder"
            return second_line,second_value,False

def check_first(first_value):
    first_line = ""
    if first_value == "Full path to the folder where the fastq files (gzipped or not) are located. Pair 1 or single-end files." or first_value == "":
        first_line = "You must put a valid path to the folder where the fastq files of the first-pair (or single-end) are located."
        return first_line,False
    else:
        return first_line,True


generateNewRun = html.Div(
    children=[
        html.H3('Required Arguments', className='required-subsection'),
        html.Div(className='box-container',
        children = [
            html.Div(
                children=[html.Label('First-pair folder',className="idminer-label"),
                    dcc.Textarea(
                        id="first-textarea",
                        className='area-text',
                        value= 'Full path to the folder where the fastq files (gzipped or not) are located. Pair 1 or single-end files.',
                        title="Full path of First-Pair Folder. i.e. /home/Maria/Experiment/FirstPair/"
            )]),
            html.Div(
                children=[html.Label('Second-pair folder',className="idminer-label"),
                    dcc.Textarea(
                        id="second-textarea",
                        className='area-text',
                        value= 'Full path to the folder where the fastq files (gzipped or not) are located. Second Pair or leave empty if single-end',
                        title="Full path of Second-Pair Folder. i.e. /home/Maria/Experiment/SecondPair/"
            )])
        ]),
        html.Div(className='box-container',
        children=[
            html.Div(
                children=[html.Label('Full path to genome file',className="idminer-label"),
                    dcc.Textarea(
                        id="genome-textarea",
                        className='area-text',
                        value= 'Fasta format',
                        title="Full path where the genome file is located. Must be in fasta format. i.e. /home/Maria/Genome/genome.fasta"
            )]),
            html.Div(
            children=[html.Label('Full path to annotation file',    className="idminer-label"),
                dcc.Textarea(
                    id="annotation-textarea",
                    className='area-text',
                    value= 'GFF format',
                    title="Full path where the annotation file is located. Must be in gff format. i.e. /home/Maria/Annotation/annotation.gff"
            )])
        ]),
        html.Div(className='box-container',
        children=[
            html.Div(
            children=[
                html.Label('Type of experiment',className="idminer-label"),
                dcc.Dropdown(
                    id='experiment-dropdown',
                    options=[
                        {'label': 'Single-end', 'value': "Single-end"},
                        {'label': 'Paired-end', 'value': 'Paired-end'},
                    ],
                    value="Paired-end",
                    placeholder='Type of experiment',
                    className="freq-dropdown"
                )]),
            html.Div(
            children=[
                html.Label('Organism',title="Define the spliced-ledear sequence used in the program.",className="idminer-label"),
                dcc.Dropdown(
                    id='organism-dropdown',
                    options=[
                        {'label': 'T. cruzi', 'value': 'T. cruzi'},
                        {'label': 'T. brucei', 'value': 'T. brucei'},
                        {'label': 'L. major', 'value': 'L. major'},
                        {'label': 'Other', 'value': 'Other'}
                    ],
                    value="T. cruzi",
                    placeholder='Select level...',
                    className="freq-dropdown"
            )]),
            html.Div(
            children=[
                html.Label('Spliced-leader sequence',className="idminer-label"),
                dcc.Textarea(
                    id='sl-textarea',
                    title='If OTHER is selected, write SL sequence of interest.',
                    className='area-text',
                    value=''
            )]),
            html.Div(
            children=[
                html.Label('Basename',className="idminer-label"),
                dcc.Textarea(
                    id='basename-textarea',
                    title='Basename of the output files.',
                    className='area-text',
                    value='UTRme-Run'
            )]),
        ]),
        
        html.Br(),
        html.H3('Optional Arguments', className='required-subsection'),
        html.Div(className='box-container',
        children=[
            html.Div(
            children=[
                html.Label('Feature type',title="Feature type (3rd column in GFF file) to be used, all features of other type are ignored.",className="idminer-label"),
                dcc.Dropdown(
                    id='feature-dropdown',
                    options=[
                        {'label': 'Gene', 'value': 'gene'},
                        {'label': 'CDS', 'value': 'CDS'},
                        {'label': 'Polypeptide', 'value': 'polypeptide'},
                        {'label': 'mRNA', 'value': 'mRNA'}
                    ],
                    value="CDS",
                    placeholder='Select feature...',
                    className="freq-dropdown"
                )]),
            html.Div(
            children=[
                html.Label('Min. overlap length',title="Select overlap length",className="idminer-label"),
                dcc.Dropdown(
                    id='overlap-dropdown',
                    options=[
                        {'label': 3, 'value': 3},
                        {'label': 4, 'value': 4},
                        {'label': 5, 'value': 5},
                        {'label': 6, 'value': 6},
                        {'label': 7, 'value': 7},
                        {'label': 8, 'value': 8},
                        {'label': 9, 'value': 9},
                        {'label': 10, 'value': 10},
                    ],
                    value=5,
                    placeholder='Select overlap length...',
                    className="freq-dropdown"
            )]),
            html.Div(
            children=[
                html.Label('Error probability',title="Cutadapt option: All searches for secondary regions are error tolerant",className="idminer-label"),
                dcc.Dropdown(
                    id='error-dropdown',
                    options=[
                        {'label': 0.05, 'value': 0.05},
                        {'label': 0.03, 'value': 0.03},
                        {'label': 0.02, 'value': 0.02},
                        {'label': 0.01, 'value': 0.01},
                        {'label': 0.005, 'value': 0.005}
                        ],
                    value=0.05,
                    placeholder='Select overlap length...',
                    className="freq-dropdown"
                    )]),
            html.Div(
            children=[
                html.Label("5'UTR length",title="Max. 5'UTRs Length",className="idminer-label"),
                dcc.Dropdown(
                    id='5len-dropdown',
                    options=[
                        {'label': 500, 'value': 500},
                        {'label': 1000, 'value': 1000},
                        {'label': 2000, 'value': 2000},
                        {'label': 3000, 'value': 3000},
                        {'label': 5000, 'value': 5000},
                        {'label': 10000, 'value': 10000},
                        {'label': "no filter", 'value': "no filter"}
                    ],
                    value=1000,
                    placeholder='Select 5 UTR length...',
                    className="freq-dropdown"
            )]),
            html.Div(
            children=[
                html.Label("3'UTR length",title="Max. 3'UTRs Length",className="idminer-label"),
                dcc.Dropdown(
                    id='3len-dropdown',
                    options=[
                        {'label': 500, 'value': 500},
                        {'label': 1000, 'value': 1000},
                        {'label': 2000, 'value': 2000},
                        {'label': 3000, 'value': 3000},
                        {'label': 5000, 'value': 5000},
                        {'label': 10000, 'value': 10000},
                        {'label': "no filter", 'value': "no filter"}
                        ],
                    value=3000,
                    placeholder='Select 3 UTR length...',
                    className="freq-dropdown"
                )]),
        html.Div(
            children=[
                html.Label("Max. ORF Length (aa)",title="Max. ORF length (aa) in UTR",className="idminer-label"),
                dcc.Dropdown(
                    id='orf-dropdown',
                    options=[
                        {'label': 30, 'value': 30},
                        {'label': 50, 'value': 50},
                        {'label': 100, 'value': 100},
                        {'label': 200, 'value': 200},
                        {'label': 300, 'value': 300},
                        {'label': "no filter", 'value': "no filter"}
                        ],
                    value=200,
                    placeholder='Remove MultiMapping Reads',
                    className="freq-dropdown"
                )]),
        html.Div(
        children=[
            html.Label('Id attribute',className="idminer-label"),
            dcc.Textarea(
                id='id-textarea',
                title='GFF attribute to be used as feature ID. i.e. ID=TcCLB.506779.120',
                className='area-text',
                value='ID='
        )]),
        html.Div(
        children=[
            html.Label('Adapter',className="idminer-label"),
            dcc.Textarea(
                id='adapter-textarea',
                title='Adapter sequences to filter out. Leave empty if none.',
                className='area-text',
                value='AGATCGGAAGAGC'
        )]),
        html.Div(
        children=[
            html.Label('Threads',className="idminer-label"),
            dcc.Textarea(
                id='core-textarea',
                title='Number of parallel search cores.',
                className='area-text',
                value=4
        )]),
        html.Div(
            children=[
                html.Label("MultiMapping",title="Remove MultiMapping Reads",className="idminer-label"),
                dcc.Dropdown(
                    id='multimapping-dropdown',
                    options=[
                        {'label': "YES", 'value': "YES"},
                        {'label': "NO", 'value': "NO"},
                        ],
                    value="YES",
                    placeholder='Remove MultiMapping Reads',
                    className="freq-dropdown"
        )]),
        html.Div(
            children=[
                html.Label("Perform analysis in 3' UTR",title="Perform analysis in 3 'UTR",className="idminer-label"),
                dcc.Dropdown(
                    id='3UTR-dropdown',
                    options=[
                        {'label': "YES", 'value': "YES"},
                        {'label': "NO", 'value': "NO"},
                        ],
                    value="YES",
                    placeholder="Perform analysis in 3' UTR",
                    className="freq-dropdown"
        )]),
        html.Div(
            children=[
                html.Label("Perform analysis in 5' UTR",title="Perform analysis in 5' UTR",className="idminer-label"),
                dcc.Dropdown(
                    id='5UTR-dropdown',
                    options=[
                        {'label': "YES", 'value': "YES"},
                        {'label': "NO", 'value': "NO"},
                        ],
                    value="YES",
                    placeholder="Perform analysis in 5' UTR",
                    className="freq-dropdown"
        )]),
        html.Div(
            children=[
                html.Label("Report UTR's with negative score",title="Show bad predicitons",className="idminer-label"),
                dcc.Dropdown(
                    id='score-dropdown',
                    options=[
                        {'label': "YES", 'value': "YES"},
                        {'label': "NO", 'value': "NO"},
                        ],
                    value="NO",
                    placeholder="Report UTR's with negative score",
                    className="freq-dropdown"
        )]),
        html.Div(
            children=[
                html.Label("Report UTR's with N's",title="Keep UTR with N in their sequence",className="idminer-label"),
                dcc.Dropdown(
                    id='N-dropdown',
                    options=[
                        {'label': "YES", 'value': "YES"},
                        {'label': "NO", 'value': "NO"},
                        ],
                    value="NO",
                    placeholder='Keep UTR with N in their sequence',
                    className="freq-dropdown"
        )]),
        html.Div(
            children=[
                html.Label("Excel",title="Write output as Excel file",className="idminer-label"),
                dcc.Dropdown(
                    id='excel-dropdown',
                    options=[
                        {'label': "YES", 'value': "YES"},
                        {'label': "NO", 'value': "NO"},
                        ],
                    value="NO",
                    placeholder='Write output as Excel file',
                    className="freq-dropdown"
        )]),
        html.Div(
            children=[
                html.Label("Delete temporary folder",title="Delete temporary folder",className="idminer-label"),
                dcc.Dropdown(
                    id='temporary-dropdown',
                    options=[
                        {'label': "YES", 'value': "YES"},
                        {'label': "NO", 'value': "NO"},
                        ],
                    value="YES",
                    placeholder='Delete temporary folder',
                    className="freq-dropdown"
        )])
    ]),
    html.Div(className="container4",
        children = [
            html.Button(
                "Create Configuration File!",
                id="run-btn",
                title="""
                Please see the wheel spinning at the bottom of the of the page.
                When the run is complete, you would see an status message in the left corner of your screen.
                """
    )])
])

layout = html.Div(
    children=[
    header.layout_home,
    html.Div(
        id='configuration-form-container',
        children=[
        html.Div(
            id='login-container-centered',
            children=[
                generateNewRun,
                html.P(id="report-status")
            ]
    )]),
    html.Div(id="return-button",style={"margin-bottom": "40px","margin-top": "30px","margin-left": "10%",    "font-size": "xx-large","font-style": "normal"})
    ])


@app.callback(dash.dependencies.Output("return-button","children"), [dash.dependencies.Input('run-btn', 'n_clicks')],[dash.dependencies.State('temporary-dropdown','value'),dash.dependencies.State('excel-dropdown','value'),dash.dependencies.State('N-dropdown','value'),dash.dependencies.State('score-dropdown','value'),dash.dependencies.State('3UTR-dropdown','value'),dash.dependencies.State('5UTR-dropdown','value'),dash.dependencies.State('multimapping-dropdown','value'),dash.dependencies.State('core-textarea','value'),dash.dependencies.State('adapter-textarea','value'),dash.dependencies.State('id-textarea','value'),dash.dependencies.State('orf-dropdown','value'),dash.dependencies.State('3len-dropdown','value'),dash.dependencies.State('5len-dropdown','value'),dash.dependencies.State('error-dropdown','value'),dash.dependencies.State('overlap-dropdown','value'),dash.dependencies.State('feature-dropdown','value'),dash.dependencies.State('basename-textarea','value'),dash.dependencies.State('sl-textarea','value'),dash.dependencies.State('organism-dropdown','value'),dash.dependencies.State('experiment-dropdown','value'),dash.dependencies.State('annotation-textarea','value'),dash.dependencies.State('genome-textarea','value'),dash.dependencies.State('second-textarea','value'),dash.dependencies.State('first-textarea','value')])
def on_click(number_of_times_button_has_clicked,temporary_value,excel_value,N_value,score_value,utr3_value,utr5_value,mmap_value,core_value,adapter_value,id_value,orf_value,len3_value,len5_value,error_value,overlap_value,feature_value,basename_value,sl_value,organism_value,experiment_value,annotation_value,genome_value,second_value,first_value):
    if number_of_times_button_has_clicked!= None:
        if temporary_value == None:
            return html.Div("You must select YES or NO in Delete temporary folder...")
        if excel_value == None:
            return html.Div("You must select YES or NO in Excel...")
        if N_value == None:
            return html.Div("You must select YES or NO in Report UTR's with N's...")
        if score_value == None:
            return html.Div("You must select YES or NO in Report UTR's with negative score...")
        if utr3_value == None:
            return html.Div("You must select YES or NO in Perform analysis in 3' UTR...")
        if utr5_value == None:
            return html.Div("You must select YES or NO in Perform analysis in 5' UTR...")
        if mmap_value == None:
            return html.Div("You must select YES or NO in MultiMapping...")
        if not check_core(core_value):
            return html.Div("You must use a integrer in Threads...")
        if not check_adapter_value(adapter_value):
            return html.Div("You must use a valid adapter sequence...")
        if not check_id(id_value):
            return html.Div("You must put an identificator...")
        if orf_value == None:
            return html.Div("You must select YES or NO in Max. ORF Length (aa)...")
        if len3_value == None:
            return html.Div("You must select a value in 3'UTR length...")
        if len5_value == None:
            return html.Div("You must select a value in 5'UTR length...")
        if error_value == None:
            return html.Div("You must select a value in Error probability...")
        if overlap_value == None:
            return html.Div("You must select a value in Min. overlap length...")
        if feature_value == None:
            return html.Div("You must select a value in Feature Type...")
        basename_value,flag = check_basename(basename_value)
        if not flag:
            return html.Div("Basename must start with a character...")
        configuration_file = basename_value + "_configuration.txt"
        sl_line,flag = check_sl_value(sl_value,organism_value)
        if not flag:
            return html.Div(sl_line)
        if organism_value == None:
            return html.Div("You must select a value in Organism...")
        if experiment_value == None:
            return html.Div("You must select a value in Type of experiment...")
        annotation_line,flag = check_annotation(annotation_value)
        if not flag:
            return html.Div(annotation_line)
        genome_line,flag = check_genome(genome_value)
        if not flag:
            return html.Div(genome_line)
        second_line,second_value,flag = check_second(second_value,experiment_value)
        if not flag:
            return html.Div(second_line)
        first_line,flag = check_first(first_value)
        if not flag:
            return html.Div(first_line)           
        display_line = "Configuration File Created! Look in Configuration_File directory."
        temporary_value,excel_value,N_value,score_value,utr3_value,utr5_value,mmap_value,core_value,adapter_value,id_value,orf_value,len3_value,len5_value,error_value,overlap_value,feature_value,basename_value,sl_value,organism_value,experiment_value,annotation_value,genome_value,second_value,first_value = map(str,[temporary_value,excel_value,N_value,score_value,utr3_value,utr5_value,mmap_value,core_value,adapter_value,id_value,orf_value,len3_value,len5_value,error_value,overlap_value,feature_value,basename_value,sl_value,organism_value,experiment_value,annotation_value,genome_value,second_value,first_value])
        configuration_line = "temporary_value:"+temporary_value+"\n","excel_value:"+excel_value+"\n","N_value:"+N_value+"\n","score_value:"+score_value+"\n","utr3_value:"+utr3_value+"\n","utr5_value:"+utr5_value+"\n","mmap_value:"+mmap_value+"\n","core_value:"+core_value+"\n","adapter_value:"+adapter_value+"\n","id_value:"+id_value+"\n","orf_value:"+orf_value+"\n","len3_value:"+len3_value+"\n","len5_value:"+len5_value+"\n","error_value:"+error_value+"\n","overlap_value:"+overlap_value+"\n","feature_value:"+feature_value+"\n","basename_value:"+basename_value+"\n","sl_value:"+sl_value+"\n","organism_value:"+organism_value+"\n","experiment_value:"+experiment_value+"\n","annotation_value:"+annotation_value+"\n","genome_value:"+genome_value+"\n","second_value:"+second_value+"\n","first_value:"+first_value
        configuration_line = "".join(configuration_line)
        if not os.path.exists("Configuration_Files"):
                os.makedirs("Configuration_Files")
        output_path = os.getcwd() + "/Configuration_Files/"
        configuration_file = output_path + configuration_file
        with open(configuration_file,"w") as fo:
            fo.write(configuration_line)
        return html.Listing(display_line+"\n"+configuration_line)
    else:
        return html.Div("Waiting...")



