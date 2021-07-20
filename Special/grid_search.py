# -*- coding: utf-8 -*-
"""
Created on 2018-10-10 16:22:28
Last Modified on 2018-10-10 16:22:28

Find best model and its parameter

@Author: Ying Huang
"""
from collections import namedtuple
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import Normalizer
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC

from .self_sklearn_objs import Log_x_1
#from .self_sklearn_objs import SelfScaleFunc

def best_SVC_multi_class(
    X,
    y,
    para_grid={
        'svc__C': [0.001, 0.01, 0.1, 1, 10, 100],
        'svc__gamma': [0.001, 0.01, 0.1, 1, 10, 100]
    },
    cv=5,
    n_jobs=1,
    scoring='f1_micro',
):

    """Grid search to find best SVC multi-classification model, and the best parameters.
    
    Note: this pipline use 'log()' and 'Normalizer()' to scale data.

    Parameters:
        X: feature matrix.
        y: labels.
    Optional parameters:
        para_grid: a set of SVM parameters for Grid search. Default use 'C', 'gamma' in a scale 0.001-100.
        cv: cross validation generater. Default use stratified k-fold cross-validation (k=5).
        n_jobs: number of threads. Default 1.
        scoring: a score to evaluate classification. Default use 'f1_micro' score.
    """

    # Split training data(X) and labels(y)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=0
    )

    # make SVC pipline
    pip = make_pipeline(
        Log_x_1(),
        Normalizer(),
        SVC(
            kernel='linear',  # choose Gaussian kernel
            cache_size=1e5, # set cache size as 100GB
            random_state=0, # set random seed
        )
    )

    # generate grid search instance
    grid = GridSearchCV(
        pip, # SVC work pipline
        param_grid=para_grid, # parameters for grid
        cv=cv, # number of fold
        n_jobs=n_jobs, # number of threads
        # use 'f1_macro' score to evaluate parameters, which means each class have equal weight
        scoring=scoring, # a score to evaluate classification
    )

    # run grid search
    grid.fit(X_train, y_train)

    # make DataFrame for parameters and test score
    df_paras = pd.concat(
        [
            pd.DataFrame(grid.cv_results_['params']),
            pd.DataFrame(grid.cv_results_['mean_test_score'], columns=['mean_test_score'])
        ],
        axis=1
    ) 

    # collect all results
    res = namedtuple(
        'Result', 
        ['grid', 'best_params', 'best_score', 'test_score', 'df_paras']
    )

    return res(
        grid,
        grid.best_params_,
        grid.best_score_,
        grid.score(X_test, y_test),
        df_paras
    )
