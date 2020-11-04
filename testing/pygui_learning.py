from dearpygui.core import *
from dearpygui.simple import *

set_main_window_size(400, 400)

with window("Tutorial",width=400, height=400, x_pos=0, y_pos=0):

    add_checkbox("Checkbox")

    add_combo(name="test", items=["hi", "hi2"])


set_primary_window("Tutorial", True)

start_dearpygui()