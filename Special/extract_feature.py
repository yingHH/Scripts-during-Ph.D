# -*- coding: utf-8 -*-
"""
Created on 2018-10-10 16:07:13
Last Modified on 2018-10-10 16:07:13

Extract features from multi-classification result by using SVC 

@Author: Ying Huang
"""
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.preprocessing import Normalizer
from sklearn.svm import SVC

from .self_sklearn_objs import Log_x_1
from .self_sklearn_objs import Pipeline_coef


def select_by_RFECV_SVC(X, y, params, step=0.01, cv=5, score='f1_micro', n_jobs=1):

    # Make pipeline
    pip = Pipeline_coef(
        [
            ("scale_log", Log_x_1()),  # log(X + 1)
            ("Scale_norm", Normalizer()),  # Normalizer().fit_tranform(X)
            ("SVC", SVC(
                kernel='linear',  # choose Gaussian kernel
                cache_size=1e5,  # set cache size as 100GB
                random_state=0,
                **params,
            )),
        ]
    )

    # Cearte object to select features using RFECV (Recursive feature elimination with cross-validation)
    rfecv = RFECV(
        pip,
        step=step,
        cv=StratifiedKFold(5),
        scoring=score,
        n_jobs=n_jobs,
    )
    
    # fit RFECV
    rfecv.fit(X, y)

    return rfecv

