#Generated Cade by DeepSeek 8/24/2025.

"""Gui Wrapper for examples
"""

import sys
import numpy as np
import pyqtgraph as pg
from PyQt6 import QtWidgets, QtCore

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PyQtGraph Application")
        self.setGeometry(100, 100, 1200, 800)

        # Central widget and main layout
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        layout = QtWidgets.QHBoxLayout(central_widget)

        # Left panel for controls
        control_panel = QtWidgets.QWidget()
        control_layout = QtWidgets.QVBoxLayout(control_panel)
        control_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)

        # Create dropdown lists with dummy data
        self.dropdowns = {}
       
        # Dropdown 1: Measurement Type
        dropdown1 = QtWidgets.QComboBox()
        dropdown1.addItems(["Voltage", "Current", "Resistance", "Power", "Frequency"])
        dropdown1.currentTextChanged.connect(self.on_measurement_changed)
        control_layout.addWidget(QtWidgets.QLabel("Measurement Type:"))
        control_layout.addWidget(dropdown1)
        self.dropdowns["measurement"] = dropdown1

        # Dropdown 2: Range Selection
        dropdown2 = QtWidgets.QComboBox()
        dropdown2.addItems(["Auto", "10 mV", "100 mV", "1 V", "10 V", "100 V"])
        dropdown2.currentTextChanged.connect(self.on_range_changed)
        control_layout.addWidget(QtWidgets.QLabel("Range:"))
        control_layout.addWidget(dropdown2)
        self.dropdowns["range"] = dropdown2

        # Dropdown 3: Sample Rate
        dropdown3 = QtWidgets.QComboBox()
        dropdown3.addItems(["1 Hz", "10 Hz", "100 Hz", "1 kHz", "10 kHz", "100 kHz"])
        dropdown3.currentTextChanged.connect(self.on_sample_rate_changed)
        control_layout.addWidget(QtWidgets.QLabel("Sample Rate:"))
        control_layout.addWidget(dropdown3)
        self.dropdowns["sample_rate"] = dropdown3

        # Create buttons
        self.buttons = {}
        button_config = [
            ("Start", self.on_start_clicked),
            ("Stop", self.on_stop_clicked),
            ("Analyze", self.on_analyze_clicked),
            ("Export", self.on_export_clicked)
        ]

        for name, callback in button_config:
            btn = QtWidgets.QPushButton(name)
            btn.clicked.connect(callback)
            self.buttons[name] = btn
            control_layout.addWidget(btn)

        # Add stretch to push content to top
        control_layout.addStretch()

        # Right panel for plot
        plot_widget = pg.GraphicsLayoutWidget()
        self.plot = plot_widget.addPlot(title="Data Plot")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.plot_curve = self.plot.plot(pen=pg.mkPen('b', width=2))

        # Sample data
        self.x_data = np.linspace(0, 10, 100)
        self.y_data = np.sin(self.x_data)
        self.plot_curve.setData(self.x_data, self.y_data)

        # Add widgets to main layout
        layout.addWidget(control_panel, 1)
        layout.addWidget(plot_widget, 4)

    # Dropdown callback functions
    def on_measurement_changed(self, text):
        print(f"Measurement type changed to: {text}")
        # Implement your measurement type change logic here

    def on_range_changed(self, text):
        print(f"Range changed to: {text}")
        # Implement your range change logic here

    def on_sample_rate_changed(self, text):
        print(f"Sample rate changed to: {text}")
        # Implement your sample rate change logic here

    # Button callback functions
    def on_start_clicked(self):
        print("Start button clicked")
        current_measurement = self.dropdowns["measurement"].currentText()
        current_range = self.dropdowns["range"].currentText()
        current_sample_rate = self.dropdowns["sample_rate"].currentText()
        print(f"Starting with: {current_measurement}, {current_range}, {current_sample_rate}")
        # Implement data acquisition/plotting start

    def on_stop_clicked(self):
        print("Stop button clicked")
        # Implement stopping procedure

    def on_analyze_clicked(self):
        print("Analyze button clicked")
        # Implement data analysis
        processed_data = self.y_data * 2  # Example transformation
        self.plot_curve.setData(self.x_data, processed_data)

    def on_export_clicked(self):
        print("Export button clicked")
        options = QtWidgets.QFileDialog.Option(0)
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save Data", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if filename:
            print(f"Exporting to {filename}")
            # Implement your export logic here

def main():
    app = QtWidgets.QApplication(sys.argv)
   
    # Optional: Set application style
    app.setStyle('Fusion')
   
    window = MainWindow()
    window.show()
   
    sys.exit(app.exec())

if __name__ == '__main__':
    main()