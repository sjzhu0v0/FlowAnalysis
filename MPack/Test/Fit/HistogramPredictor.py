import numpy as np
import matplotlib.pyplot as plt
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from typing import Callable

class HistogramPredictor:
    def __init__(self):
        self.model = None
        self.input_shape = None

    def train(self, func: Callable, param_range, n_histograms, n_samples=1000, n_bins=30, test_size=0.2, epochs=100, batch_size=32):
        X_train, X_test, y_train, y_test = generate_data(func, param_range, n_histograms, n_samples=n_samples, n_bins=n_bins)
        self.input_shape = (X_train.shape[1],)
        self.model = build_model(self.input_shape)
        self.model.fit(X_train, y_train, validation_data=(X_test, y_test), epochs=epochs, batch_size=batch_size)

    def predict(self, histograms):
        if self.model is None:
            raise Exception("Model not trained yet. Call `train` method before predicting.")
        return self.model.predict(histograms)

    def get_trained_model(self):
        if self.model is None:
            raise Exception("Model not trained yet. Call `train` method before getting the model.")
        return self.model

