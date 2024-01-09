from ema import EyeMovementData
from ema import Model
from ema import Htmlreport

data = EyeMovementData()

model = Model(data, k=3, dynamic_model=True, dynamic_model_frequency=5)
model.iterate_em(500)

html_report = Htmlreport(model, number_of_text_restorations_to_display=0)
html_report.make_html(open_in_web_browser=True)
