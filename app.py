import dash
from flask import Flask, send_from_directory



# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:

server = Flask(__name__)


app = dash.Dash(__name__,server=server)
app.config.supress_callback_exceptions = True


# Add bootstrap css

external_css = [
    'https://codepen.io/chriddyp/pen/bWLwgP.css',
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
    "//fonts.googleapis.com/css?family=Raleway:400,300,600",
    "https://codepen.io/bcd/pen/KQrXdb.css",
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css",
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css"
    ]

for css in external_css:
    app.css.append_css({"external_url": css})
