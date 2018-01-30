EMEWS project template
-----------------------

You have just created an EMEWS project.

This project is compatible with swift-t v. 1.3+. Earlier
versions will NOT work.

The project consists of the following directories:

```
  data/
  ext/
  etc/
  scripts/
  swift/
  model/
  README.md
```
The directories are intended to contain the following:

 * `data` - model input etc. data
 * `etc` - additional code used by EMEWS
 * `scripts` - any necessary scripts (e.g. scripts to launch a model), excluding
    scripts used to run the workflow.
 * `swift` - swift code
 * `model` - should contain the model itself

The compiled PhysiCell model should be copied into the model directory and be named
cancer-immune-EMEWS.
