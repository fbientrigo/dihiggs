import unittest
import logging
import time

# chequear los improts y modificar como ocurre el testing
from parallel_utils import OracleExecutor

# Configure root logger to output debug messages to console
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s [%(levelname)s] %(message)s')

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
