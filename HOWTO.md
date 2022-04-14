OXSX/OXO is a signal extraction frmaework built for the SNO+ experiment.
It is built to be useful for general analyses within particle physics.

This document will try and explain the critical components of OXO, so you can
build your analysis code efficiently.

- [1. Objects for storing data](#1-objects-for-storing-data)
  - [1.1. Event](#11-event)
  - [1.2. ObsSet](#12-obsset)
  - [1.3. DataSet](#13-dataset)
    - [1.3.1. OXSXDataSet](#131-oxsxdataset)
    - [1.3.2. ROOTNtuple](#132-rootntuple)
    - [1.3.3. LazyOXSXDataSet](#133-lazyoxsxdataset)

## 1. Objects for storing data
### 1.1. Event

The `Event` class repesents the observations associated with a single event.
For example, you might have observed an event with properties:
 - Reconstructed energy = 2.54 MeV
 - Reconstructed radius = 3.17 m
 - External classifier value = 0.6

The event class would then store this information as a vector of information
(2.54, 3.17, 0.6), along with a vector of the same size indicating the names
for the associated observables, `("energy", "radius", "ext_classifier")`, say.

In general, an `Event` object stores the data for any set of observables given
to it, as long as the observables are doubles.

Creating & handling `Event` objects is easy (though you'll probably handle individual
events on their own very rarely in your analysis, I suspect). For our example,
we could create an `Event` object like so:
```
const std::vector<double> event_data {2.54, 3.17, 0.6};
const std::vector<std::string> observable_names {"energy", "radius", "ext_classifier"};
Event ev(event_data);
ev.SetObservableNames(&observable_names);
```
Once an `Event` object has been created, you can access the information inside
it, via `GetDatum()` for a single observable value, or `GetData()` and
`GetDataMap()` for all the observables at once. `GetDatum()` can be used 
either
with the observable name, or the index of the observables vector.

### 1.2. ObsSet

`ObsSet` is like `Event`, just without any of the actual data! It's just a
store of the observable names on their own.

Given an `Event` `ev`, one can get the data for a particular set of 
observables
defined in an `ObsSet` `obs` via:
```
const std::vector<double> obs_data = ev.(obs);
```

### 1.3. DataSet

`DataSet` is just an abstract base class for some much more useful classes
we'll talk about in a moment. The point of any such object is to store a
bunch of data. The most important method here is `GetEntry(i)`, where you
get the ith event in the dataset.

Here are the derived classes of `DataSet`, which can actually be used:

#### 1.3.1. OXSXDataSet

The most obvious kind of `DataSet` object in OXO: it holds a vector of `Events` (and the assocaited observable names in another vector).
In addition to `GetEntry(i)`, one can also add an event to the dataset by
`AddEntry()`. A nice feature is concatenating `OXSXDataSet` objects together
with the `+` operator overload. One can reserve memory for n events by 
`Reserve(n)`.

#### 1.3.2. ROOTNtuple

A `DataSet` object derived from a ROOT `TNtuple` object that is stored in a
file. When the object is created, the `TNtuple` object is loaded. No adding
events afterwards is then possible.

#### 1.3.3. LazyOXSXDataSet

A fancy version of `OXSXDataSet`; so fancy I'm not sure anyone uses this.
You create the object with a filename pointing to some data that is stored
within a `.H5` file (more on this in the IO section). Getting data happens
as usual with `GetEntry(i)`. There's some shenanigans with a static list of
these objects stored, leading to the functionality where if you have multiple
of these objects at once I'm pretty sure only one can have data loaded at a time.

# Histograms, etc.
Histograms are really useful in Particle Physics analyses! But to describe & use
them in the most general terms, we'll first have to create some other objects 
that will help us.
## BinAxis

`BinAxis` is just a class that stores the bin edges for some variable. It can be
set up in two ways: either with equal bin widths, or unequal widths. In the
latter case, the user must specify the bin edges. The axis must come along with
a name (e.g. "energy"), and can also have an optional name to be used for pretty LaTeX formatting ("E_{reco}", for example).

You can do most of the usual things you would expect for a binned axis with this
object. You will hopefully not have to interact too much with this class, other
than to set up axes for any histogramming you do.

## AxisCollection

As the name suggests, this class holds a number of `BinAxis` objects. To
actually create one with the binning you want, first create an empty 
`AxisCollection` object, and then use either `AddAxis()` or `AddAxes()` as
appropriate.

Technically, you can call a bunch of methods with this class. Most of them you'll
never need to call yourself. Like `BinAxis` and `ObsSet`, `AxisCollection` objects mostly have to be set up at the beginning of your code, and then never really
have to be touched.

## Histogram

This class stores an N-dimensional histogram, the binning of which is defined
by the `AxisCollection` object you give it.

Like a ROOT `TH1`, the `Fill()` method is used to add data to the histogram 
(you can add values one-by-one, or give a whole vector of them at once), which
can be weighted. You can also modify the values of the bins directly through
`SetBinContent()`.

This class has a whole bunch of nice things you can do with it! You can normalise,
scale, and get the Integral of the histogram; 