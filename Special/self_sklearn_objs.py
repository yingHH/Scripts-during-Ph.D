# -*- coding: utf-8 -*-
"""
Created on 2018-10-12 19:46:13
Last Modified on 2018-10-12 19:46:13

Create self-sklearn objects by inheriting base class 'sklearn.base' 

@Author: Ying Huang
"""
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline


# !!! Have some problem in 'scale_func'
class SelfScaleFunc(BaseEstimator, TransformerMixin):
    """
    Scaling trianing data by using self-scaling function.

    Parameters:
        scale_func: self-funtion to scale training data.

    Example:
        # Scaling training data (X) by log(x + 1).
        import numpy as np


        def log_x_1(x):
            return np.log(x + 1)

        scale_log = SelfScaleFunc(log_x_1)
        scale_log.fit(X)
        scale_log.transform()
    """

    def __init__(self, scale_func):
        self.scale_func = scale_func

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):

        try:
            if y is None:
                print(X)
                return self.scale_func(X)
            elif y is not None:
                return self.scale_func(X, y)
        except:
            msg = """
            ERR: 'scale_func' is not working, please check your scale function.
            Note: it is better to use pakages 'Numpy' or 'Pandas' to process data."""
            raise Exception(msg)


class Log_x_1(BaseEstimator, TransformerMixin):
    """
    Scaling training data (X) by log(x + 1).

    Parameters:
        X: training data.

    Example:

        scale_log = Log_x_1()
        scale_log.fit(X)
        scale_log.transform(X)
    """

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return np.log(X + 1)

    def fit_transform(self, X, y=None, **fit_params):
        return np.log(X + 1)


class Pipeline_coef(Pipeline):

    @property
    def coef_(self):
        # Get final step name. Usually it is a name of estimator
        final_step_name = self.steps[-1][0]
        # Get the coef_ feature of the estimator(final step), and return it
        try:
            return self.named_steps[final_step_name].coef_
        except:
            msg = """
            Find no 'coef_' feature. Please check below:
            * The final step should be a estimator that contains 'coef_' feature.
            * 'coef_' can be called after running Pipeline_coef(...).fit(X, y).
            """
            raise Exception(msg)
