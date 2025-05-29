# Advanced appliactions
This part of the library depends on source code at the base, creating interconnected modular abstraction of code, capable of clearly be used in dynamic environments.


## Info Model
Allows to input both the param space and the code to be run, accepting a function as input and working as a functional that 

```
def associated_model:

    def __init__.py(self, function):
        """
        the function must:
            accept a paramer space, and **kwargs in order to use a normal executtion or parallel execution.
        """
        ...
        # a try or test its done in order to test the function:
        ...
        ## callables are constructed
        if non_parable:
            ....

        # the functons get saved in the model
        self.model = function

    def infer_time_space(self, search_space : np.array ):
        """
        answer data about time and space when calcullating a given associated parameter space
        """
        ...
```


## depends
- `utils.py`:
    - ``
    - `estimate_exploration_time`



## Threading
Allows to execute parallel code.

### depends:
- oracle.py
    - `def run_oracle_batch(param_list, nthreads, executable_path=EXECUTABLE_PATH, debug=False):`
    - `def safe_run_oracle(params):` 

### Use case:

...

## Plotting
Accepts different plot styles inside a class object. Capable of having the data loaded on non RAM sectors temporally and then used in liner $\mathcal O$ scale

### depends:
- 


