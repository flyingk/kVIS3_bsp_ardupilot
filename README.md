## kVIS3 Board Support Package (BSP) for ArduPilot Data

A support package for [kVIS3 data visualisation](https://github.com/flyingk/kVIS3).

### Features:
- import an ArduPilot dataset (*.bin)
- plot flight path on google map (API key required)
- custom plot definitions

### Installation:
#### Getting the Code
The ArduPilot BSP can be cloned using git into your preferred directory, or downloaded as a zip from the github page.
```
git clone https://github.com/flyingk/kVIS3_bsp_ardupilot.git --recursive
```
The import relies on a submodule to import the BIN file into MATLAB, so if you forget the `--recursive` tag, you can pull the submodules from within the BSP but running
```
git submodule init
git submodule update
```

#### Loading the BSP into kVIS
The BSP needs to be set in kVIS.  This can be accomplished by running kVIS3, navigating to `Edit > Change BSP`, and selecting the folder for the BSP.

### Contributing
Issues, pull requests, and contributions are always welcome!

Thanks for watching, starring and sharing!
