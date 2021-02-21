import plotly.offline
import plotly.graph_objs as go
import plotly.figure_factory as ff
import numpy as np
from sklearn import metrics

def create_probability_histogram(correct_probabilities, incorrect_probabilities, save_to_file=False):
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=correct_probabilities, name="Correct predictions", nbinsx=100, marker_color='#009988',
                               hovertemplate = '<b>Range</b>: %{x}<br>'+ '<b>Predictions</b>: %{y}'))
    fig.add_trace(go.Histogram(x=incorrect_probabilities, name="Incorrect predictions", nbinsx=100, marker_color='#CC3311',
                               hovertemplate = '<b>Range</b>: %{x}<br>'+ '<b>Predictions</b>: %{y}'))

    fig.update_layout(barmode='overlay', xaxis={"title": "Prediction certainty"}, yaxis={"title": "Prediction count"},
                      margin=go.layout.Margin(t=35, l=50, r=200, autoexpand=False), template='plotly_white')
    fig.update_traces(opacity=0.75)
    fig.update_xaxes(nticks=11, range=[0, 1])
    fig.update_layout()

    if save_to_file:
        plotly.offline.plot(fig, filename='histogram.html', auto_open=False)

    fig.show()
    
    
def plot_calibration(correct_probabilities, incorrect_probabilities, save_to_file=False):
    c1, b1 = np.histogram(correct_probabilities, bins=100)
    c2, b2 = np.histogram(incorrect_probabilities, bins=100)

    fig = go.Figure()
    fig.add_trace(go.Scatter(name="CMB classifier", x=b1[:].round(4)*100, y=c1/(c1+c2).astype(float).round(4)*100,
                             hovertemplate = '<b>Certainty range</b>: %{x:.2f}%<br>'+ 
                             '<b>Prediction accuracy</b>: %{y:.2f}%', marker_color='#1EC283'))
    fig.add_trace(go.Scatter(name="Perfect calibration", x=b1[:].round(4)*100, y=b1[:].round(4)*100, marker_color='#999999',
                             hovertemplate = '<b>Certainty range</b>: %{x:.2f}%<br>'+ 
                             '<b>Prediction accuracy</b>: %{y:.2f}%'))

    fig.update_layout(xaxis={"title": "Prediction certainty [%]"}, yaxis={"title": "Prediction accuracy [%]"},
                      margin=go.layout.Margin(t=50, l=50, r=200, autoexpand=False), template='plotly_white')
    fig.update_yaxes(nticks=11, range=[0, 100])
    fig.update_xaxes(nticks=11, range=[0, 100])

    if save_to_file:
        plotly.offline.plot(fig, filename='calibration.html', auto_open=False)

    fig.show()


def plot_interactive_confusion_matrix(predictions, sort_by_popularity=False, save_to_file=False):
    if sort_by_popularity:
        class_labels = list(predictions.y_true.value_counts().index)
    else:
        class_labels = sorted(list(predictions.y_true.unique()))
        
    cm = metrics.confusion_matrix(predictions.y_true, predictions.y_pred, labels=class_labels)
    
    cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    data = [
        go.Heatmap(
            x=class_labels,
            y=class_labels[::-1],
            z=cm_normalized[::-1, :].round(4)*100,
            colorscale=[[0.0, "rgb(255, 255, 255)"], [1.0, "rgb(0, 0,0)"]],
            hovertemplate =
            '<b>True ligand</b>: %{y}<br>'+
            '<b>Predicted</b>: %{x}<br>'+
            '<b>Proportion</b>: %{z:.2f}%<extra></extra>'
        )
    ]

    layout = go.Layout(
        xaxis={"title": "Predicted label"},
        yaxis={"title": "True label"},
        width=1000,
        height=775,
        autosize=False,
        margin=go.layout.Margin(
            t=15,
            l=200,
            r=200,
            autoexpand=False
        ),
    )

    fig = go.Figure(data=data, layout=layout)

    if save_to_file:
        plotly.offline.plot(fig, filename='cf.html', auto_open=False)

    fig.show()
