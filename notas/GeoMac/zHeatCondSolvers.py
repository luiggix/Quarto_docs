from cond_calor_solvers import solucion
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

metodos = ['LU (linalg)', 'Jacobi', 'Gauss-Seidel', 'SOR', 'Steepest', 'CGM']

w_metodo = widgets.Dropdown(
    options=metodos, 
    description="Método",
    layout=widgets.Layout(height='auto', width='auto')
)

w_N = widgets.IntSlider(
    min=5, max=18, step=1, value=3,
    description='N',
    layout=widgets.Layout(width='300px')
)

button = widgets.Button(
    description="", icon='play',
    layout=widgets.Layout(width='50px')    
)

ui = widgets.VBox([widgets.HBox([w_metodo, w_N, button])])

ui.layout = widgets.Layout(border='solid 1px black')
ui.layout.width = '500px'

output = widgets.Output()

display(ui, output)

def on_button_clicked(b):
    output.clear_output(wait=True)
    with output:
        solucion(B = 100, T =0, L=34, R = 7, metodo = w_metodo.value, N = w_N.value)

button.on_click(on_button_clicked)