#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `vibrational_analysis_step` package."""

import pytest  # noqa: F401
import vibrational_analysis_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = vibrational_analysis_step.VibrationalAnalysis()
    assert (
        str(type(result))
        == "<class 'vibrational_analysis_step.vibrational_analysis.VibrationalAnalysis'>"  # noqa: E501
    )
