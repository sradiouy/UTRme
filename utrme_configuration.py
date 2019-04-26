import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
import structure

#####


# loads different apps on different urls

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/':
        return structure.layout
    else:
        return 'The page does not exist: 404 Error'


if __name__ == '__main__':
    app.run_server(debug=True)