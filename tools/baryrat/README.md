## Brief
An interface for using the python library [Baryrat](https://github.com/c-f-h/baryrat).

## Install

**Step 1.** We recommend using a virtual environment. Here, we use a conda virtual environment named `env_baryrat`:

```
$ conda create -n env_baryrat python=3.8
$ conda activate env_baryrat
```

**Step 2.** The `requirements.txt` provides an easy installation using `pip`:

```
(env_baryrat) $ cd $STATFEMROOT/tools/baryrat/baryrat-x.x.x
(env_baryrat) $ pip install -r requirements.txt
```

Use `conda list` to check if baryrat is installed:

```
(env_baryrat) $ conda list | grep baryrat
```

## Call from C++

- Refer to the examples in `$STATFEMROOT/tools/baryrat/examples`. Compile and run the C++ program as usual. 

- Run in terminal: remember to activate the virtual environment **and** have the environment variable `STATFEMROOT` set.

- Run in Eclipse: set the environment variables `PYTHONHOME` and `STATFEMROOT` in your `Run` configuration setting. First, choose `Edit` > `Environment` > `Add`. After setting the environment variables, tick the option `Append environment to native environment`.
