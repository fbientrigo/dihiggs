# lib/advanced/threading.py
"""
Threading module: provides a universal Executor that can run a given model function
sequentially or in parallel using the existing `run_oracle` and `run_oracle_batch` interfaces.
"""
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable, List, Any, Dict, Optional

# chequear los imports
from ..oracle import run_oracle, run_oracle_batch, EXECUTABLE_PATH


# =========== Testing ================

import unittest
import logging
import time


# # Configure root logger to output debug messages to console
# logging.basicConfig(level=logging.DEBUG,
#                     format='%(asctime)s [%(levelname)s] %(message)s')


class OracleExecutor:
    """
    Universal executor for Oracle-based models.
    Accepts a model function that takes a parameter list and returns a dict.
    Provides sequential and parallel execution methods with standardized output.

    Attributes:
        model_fn: Callable[[List[float]], Dict]  # function for single invocation
        batch_fn: Callable[[List[List[float]], int], List[Dict]] # batch invocation
        nthreads: int  # threads to use in batch_fn
    """
    def __init__(self,
                 model_fn: Callable[[List[float]], Dict] = run_oracle,
                 batch_fn: Callable[[List[List[float]], int], List[Dict]] = run_oracle_batch,
                 nthreads: int = 4):
        self.model_fn = model_fn
        self.batch_fn = batch_fn
        self.nthreads = nthreads

    def run(self, params: List[float]) -> Dict:
        """
        Run single parameter set sequentially.
        Returns a standardized dict with output or error.
        """
        try:
            result = self.model_fn(params)
            return self._standardize(result)
        except Exception as e:
            return self._error_dict(params, str(e))

    def run_batch(self, param_list: List[List[float]]) -> List[Dict]:
        """
        Run multiple parameter sets in batch using the C++ parallel backend.
        Each param set must be length 7.
        """
        try:
            results = self.batch_fn(param_list, self.nthreads)
        except Exception as e:
            # on batch error, fallback to sequential
            results = [self.run(p) for p in param_list]
        return [self._standardize(r) for r in results]

    def map(self, param_list: List[List[float]], use_threads: bool = False) -> List[Dict]:
        """
        Generic map: if use_threads=False uses batch backend, if True uses ThreadPoolExecutor
        for individual threads.
        """
        if use_threads:
            with ThreadPoolExecutor(max_workers=self.nthreads) as exe:
                futures = [exe.submit(self.run, p) for p in param_list]
                return [f.result() for f in futures]
        else:
            return self.run_batch(param_list)

    def _standardize(self, output: Dict) -> Dict:
        """
        Ensure all expected keys present, even on error.
        """
        # expected_keys similar to safe_run_oracle
        expected = {
            'positivity_ok': None,
            'unitarity_ok': None,
            'perturbativity_ok': None,
            'w_h2_bb': None,
            'w_h2_tautau': None,
            'w_h2_uu': None,
            'w_h2_du': None,
            'w_h2_ln': None,
            'w_h2_vv': [None, None, None],
            'w_h2_gaga': None,
            'w_h2_Zga': None,
            'w_h2_gg': None,
            'w_h2_hh': None,
            'w_total_h2': None,
            'w_total_top': None,
            'branching_ratio_h2_gaga': None,
            'lambda1': None,
            'lambda2': None,
            'lambda3': None,
            'lambda4': None,
            'lambda5': None,
            'lambda6': None,
            'lambda7': None,
        }
        std = expected.copy()
        # If error in output, preserve it
        if 'error' in output:
            std['error'] = output.get('error')
        for k in expected:
            std[k] = output.get(k, expected[k])
        return std

    def _error_dict(self, params: List[float], msg: str) -> Dict:
        """
        Build an error dictionary for a failed call.
        """
        return {'error': msg, 'params': params}



class DummyModel:
    """Dummy model functions for testing."""
    @staticmethod
    def good_fn(params):
        # Simulate some work
        time.sleep(0.01)
        return {'positivity_ok': 1, 'unitarity_ok': 1, 'perturbativity_ok': 1}

    @staticmethod
    def bad_fn(params):
        # Simulate failure
        raise RuntimeError("Dummy failure")

class DummyBatch:
    """Dummy batch function for testing."""
    @staticmethod
    def good_batch(param_list, nthreads):
        time.sleep(0.02)
        return [{'lambda1': i} for i, _ in enumerate(param_list)]

    @staticmethod
    def bad_batch(param_list, nthreads):
        # Simulate batch failure
        raise RuntimeError("Batch failure")

class TestOracleExecutor(unittest.TestCase):
    def setUp(self):
        # Called before each test
        logging.info("Setting up OracleExecutor for test")
        self.inputs = [[1,2,3,4,5,6,7], [7,6,5,4,3,2,1]]

    def test_run_success(self):
        logging.info("Testing run() with good_fn")
        exec = OracleExecutor(model_fn=DummyModel.good_fn)
        for params in self.inputs:
            out = exec.run(params)
            self.assertIn('positivity_ok', out)
            self.assertEqual(out['positivity_ok'], 1)
        logging.info("test_run_success completed")

    def test_run_error(self):
        logging.info("Testing run() with bad_fn")
        exec = OracleExecutor(model_fn=DummyModel.bad_fn)
        out = exec.run(self.inputs[0])
        self.assertIn('error', out)
        self.assertIn('Dummy failure', out['error'])
        logging.info("test_run_error completed")

    def test_run_batch_success(self):
        logging.info("Testing run_batch() with good_batch")
        exec = OracleExecutor(model_fn=DummyModel.good_fn,
                              batch_fn=DummyBatch.good_batch,
                              nthreads=2)
        results = exec.run_batch(self.inputs)
        self.assertEqual(len(results), len(self.inputs))
        for i, r in enumerate(results):
            self.assertEqual(r.get('lambda1'), i)
        logging.info("test_run_batch_success completed")

    def test_run_batch_fallback(self):
        logging.info("Testing run_batch() fallback on bad_batch")
        exec = OracleExecutor(model_fn=DummyModel.good_fn,
                              batch_fn=DummyBatch.bad_batch,
                              nthreads=2)
        results = exec.run_batch(self.inputs)
        # Should fallback to sequential good_fn, so positivity_ok == 1
        for r in results:
            self.assertIn('positivity_ok', r)
            self.assertEqual(r['positivity_ok'], 1)
        logging.info("test_run_batch_fallback completed")

    def test_map_threads(self):
        logging.info("Testing map(use_threads=True)")
        exec = OracleExecutor(model_fn=DummyModel.good_fn,
                              batch_fn=DummyBatch.good_batch,
                              nthreads=3)
        results = exec.map(self.inputs, use_threads=True)
        self.assertEqual(len(results), len(self.inputs))
        for r in results:
            self.assertEqual(r['positivity_ok'], 1)
        logging.info("test_map_threads completed")

    def test_map_batch(self):
        logging.info("Testing map(use_threads=False)")
        exec = OracleExecutor(model_fn=DummyModel.good_fn,
                              batch_fn=DummyBatch.good_batch,
                              nthreads=3)
        results = exec.map(self.inputs, use_threads=False)
        self.assertEqual(len(results), len(self.inputs))
        for i, r in enumerate(results):
            self.assertEqual(r.get('lambda1'), i)
        logging.info("test_map_batch completed")


if __name__ == "__main__":
    unittest.main(verbosity=2)
