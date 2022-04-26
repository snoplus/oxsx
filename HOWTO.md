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
- [2. Histograms, etc.](#2-histograms-etc)
  - [2.1. BinAxis](#21-binaxis)
  - [2.2. AxisCollection](#22-axiscollection)
  - [2.3. Histogram](#23-histogram)
- [3. PDFs & Event Distributions](#3-pdfs--event-distributions)
  - [3.1. ParameterDict](#31-parameterdict)
  - [3.2. FitComponent](#32-fitcomponent)
    - [3.2.1. Function](#321-function)
      - [3.2.1.1. Heaviside](#3211-heaviside)
      - [3.2.1.2. PDF](#3212-pdf)
        - [3.2.1.2.1. Gaussian](#32121-gaussian)
    - [3.2.2. ConditionalPDF](#322-conditionalpdf)
      - [3.2.2.1. JumpPDF](#3221-jumppdf)
      - [3.2.2.2. VaryingCDF](#3222-varyingcdf)
  - [3.3. EventDistribution](#33-eventdistribution)
    - [3.3.1. AnalyticED](#331-analyticed)
    - [3.3.2. BinnedED](#332-binneded)
    - [3.3.3. CompositeED](#333-compositeed)
    - [3.3.4. SpectralFitDist](#334-spectralfitdist)
- [4. Random number & Event generation](#4-random-number--event-generation)
  - [4.1. DataSetGenerator](#41-datasetgenerator)
  - [4.2. BinnedEDGenerator](#42-binnededgenerator)
- [5. Data type conversions & IO](#5-data-type-conversions--io)
  - [5.1. DistTools](#51-disttools)
  - [5.2. DistFiller](#52-distfiller)
  - [5.3. IO](#53-io)
- [6. Applying cuts to data](#6-applying-cuts-to-data)
  - [6.1. Cut](#61-cut)
    - [6.1.1. BoolCut](#611-boolcut)
    - [6.1.2. LineCut](#612-linecut)
    - [6.1.3. BoxCut](#613-boxcut)
  - [6.2. CutCollection](#62-cutcollection)
  - [6.3. CutLog](#63-cutlog)
- [7. Systematics](#7-systematics)
  - [7.1. EventSystematic](#71-eventsystematic)
    - [7.1.1. EventShift](#711-eventshift)
    - [7.1.2. EventScale](#712-eventscale)
    - [7.1.3. EventConvolution](#713-eventconvolution)
    - [7.1.4. EventReconvolution](#714-eventreconvolution)
  - [7.2. Systematic](#72-systematic)
    - [7.2.1. Scale](#721-scale)
    - [7.2.2. Convolution](#722-convolution)
- [Managers](#managers)
  - [ParameterManager](#parametermanager)

# 1. Objects for storing data
## 1.1. Event

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

## 1.2. ObsSet

`ObsSet` is like `Event`, just without any of the actual data! It's just a
store of the observable names on their own.

Given an `Event` `ev`, one can get the data for a particular set of 
observables
defined in an `ObsSet` `obs` via:
```
const std::vector<double> obs_data = ev.(obs);
```

## 1.3. DataSet

`DataSet` is just an abstract base class for some much more useful classes
we'll talk about in a moment. The point of any such object is to store a
bunch of data. The most important method here is `GetEntry(i)`, where you
get the ith event in the dataset.

Here are the derived classes of `DataSet`, which can actually be used:

### 1.3.1. OXSXDataSet

The most obvious kind of `DataSet` object in OXO: it holds a vector of `Events` (and the assocaited observable names in another vector).
In addition to `GetEntry(i)`, one can also add an event to the dataset by
`AddEntry()`. A nice feature is concatenating `OXSXDataSet` objects together
with the `+` operator overload. One can reserve memory for n events by 
`Reserve(n)`.

### 1.3.2. ROOTNtuple

A `DataSet` object derived from a ROOT `TNtuple` object that is stored in a
file. When the object is created, the `TNtuple` object is loaded. No adding
events afterwards is then possible.

### 1.3.3. LazyOXSXDataSet

A fancy version of `OXSXDataSet`; so fancy I'm not sure anyone uses this.
You create the object with a filename pointing to some data that is stored
within a `.H5` file (more on this in the IO section). Getting data happens
as usual with `GetEntry(i)`. There's some shenanigans with a static list of
these objects stored, leading to the functionality where if you have multiple
of these objects at once I'm pretty sure only one can have data loaded at a time.

# 2. Histograms, etc.
Histograms are really useful in Particle Physics analyses! But to describe & use
them in the most general terms, we'll first have to create some other objects 
that will help us.
## 2.1. BinAxis

`BinAxis` is just a class that stores the bin edges for some variable. It can be
set up in two ways: either with equal bin widths, or unequal widths. In the
latter case, the user must specify the bin edges. The axis must come along with
a name (e.g. "energy"), and can also have an optional name to be used for pretty LaTeX formatting ("E_{reco}", for example).

You can do most of the usual things you would expect for a binned axis with this
object. You will hopefully not have to interact too much with this class, other
than to set up axes for any histogramming you do.

## 2.2. AxisCollection

As the name suggests, this class holds a number of `BinAxis` objects. To
actually create one with the binning you want, first create an empty 
`AxisCollection` object, and then use either `AddAxis()` or `AddAxes()` as
appropriate.

Technically, you can call a bunch of methods with this class. Most of them you'll
never need to call yourself. Like `BinAxis` and `ObsSet`, `AxisCollection` objects mostly have to be set up at the beginning of your code, and then never really
have to be touched.

## 2.3. Histogram

This class stores an N-dimensional histogram, the binning of which is defined
by the `AxisCollection` object you give it.

Like a ROOT `TH1`, the `Fill()` method is used to add data to the histogram 
(you can add values one-by-one, or give a whole vector of them at once), which
can be weighted. You can also modify the values of the bins directly through
`SetBinContent()`.

This class has a whole bunch of nice things you can do with it! You can normalise,
scale, and get the integral of the histogram; get a slice of the histogram or 
marginalise over certain axes. You can also add, multiply, or divide two histograms
together. One other useful method is `AddPadding()`: this adds a tiny non-zero
value to all bins with zero in them. This comes in handy when calculating
likelihoods, as this will prevent the code from breaking entirely.

# 3. PDFs & Event Distributions
Okay, so we have classes for storing raw event information, and one for a
general binned histograms. Great! What would be nice now are some classes related
to probability and event distributions. We're gonna have to build a number of
other classes first to get there, though!

## 3.1. ParameterDict
As the name suggests, this is a dictionary of parameters, a general enough
object that it gets used all over the place in OXO. You got a bunch of named
parameters (e.g. nuisance parameters in a fit)? Put 'em in a `ParameterDict`! 
What about the means & standard deviations for a multi-dimensional Gaussian?
You can store them as a `ParameterDict`! Need to describe some settings

This class is literally just a mapping from strings --> doubles.

## 3.2. FitComponent
This is an abstract base class for objects that have parameters associated with
them. Anything derived from this class have some set of named parameters, each
of which is a double.

Any `FitComponent`-derived object can have parameters set & gotten (via a 
`ParameterDict` if you like!), and renamed. The object also has its own name
that can be changed.

### 3.2.1. Function
OXO's capital-F `Function` class is another abstract class, derived from
`FitComponent`. As the name suggests, any such object acts like a kind of 
mathematical function: it can be *evaluated* at some set of values given to it,
simply by doing `func(vals)`, where func is the `Function` object, and `vals`
is a vector of doubles.

#### 3.2.1.1. Heaviside
We finally come to our first non-abstract class in this section: `Heaviside`.
This is the multi-dimensional version of the basic mathematical function also
called the [Heaviside step function](https://mathworld.wolfram.com/HeavisideStepFunction.html).
More precisely, it returns the product of these step functions over the 
dimensions it is defined on, with a possible position of the step for each 
dimension along whther that step goes up or down (or no step for that axis).
Evaluation of this function will always return either a 1 or a zero only.

Frankly not really used very much, but it is a non-trivial example of a function
that cannot be used as a probability density function. And on that note...
#### 3.2.1.2. PDF
Let's go another level deeper! a `PDF` is *another* abstract class, which is a
`Function` that also has the property of definite integration, via the 
`Integrate(mins, maxs)` method. Because a `PDF` has a well-defined integral,
it can also be also randomly sampled from: `Sample()`. See the section on 
[random number generation](#4-random-number-generation-rand) for how derived
classes implement this.

##### 3.2.1.2.1. Gaussian
For all of that effort, we can finally make some actual PDFs, in particular a
multi-dimensional Gaussian (given how useful they are!). Unsurprisingly, 
`Gaussian` derives from the `PDF` class.

This class has all of its derived features. You can get/set the
means and standard deviations of any of the dimensions by either `SetMean()` etc.,
or with the `FitComponent` `SetParameter()` method. For the latter, the name of
the mean & sigma parameters are given by `means_i` & `stddevs_i` respectively, 
where `i` here corresponds to the index of the dimension of interest.
(Note - the actual parameter names are handles not within the `Gaussian` class
itself, but within a helper class called `GaussianFitter`. Despite the name,
this latter class doesn't fit a Gaussian, but instead holds the fit parameters.)
You can evaluate the `Gaussian` at any point in the n-dimensional space, and 
even integrate within some (possibly high-dimensional) rectangular region.

As a bonus, you can obtain the cumulative density function by `Cdf()`.

### 3.2.2. ConditionalPDF
**TODO**

#### 3.2.2.1. JumpPDF
**TODO**
#### 3.2.2.2. VaryingCDF
**TODO**

## 3.3. EventDistribution
`EventDistribution` is an abstract base class, whose derived classes describe
the expected distributions of events under a set of observables. A typical
physics analysis might involve comparison/fitting of an `EventDistribution`
expected from theory/MC, to an observed `DataSet` object.

The virtual methods any `EventDistribution` object can perform include 
calculating the probability that a given `Event` object arises from the event
distribution; the total integral of the event distribution; and a `Normalise()`
method which will scale the event distribution so that the integral is 1.
Finally, the event distribution has a name that can be set/gotten.

Note: it is this class that Particle Physicists typically refer to as "PDFs".
We're often using that terminology fast and loose, given that true PDFs must
have a total normalisation of 1 by definition. Because of this, actual
mathematical PDFs are defined within the `PDF` class, separate from this class.
### 3.3.1. AnalyticED
One possible event distribution is a purely analytic one, givne by `AnalyticED`.
This class holds a `PDF` that holds the underlying shape of the distribution,
along with a constant normalisation. It also has an `ObsSet` object, so that
the observable quantities of interest can be defined distinct to the (possibly)
larger set of variables within the PDF. On top of this, `AnalyticED` inherits 
from `FitComponent`, with the parameters being defined within the `PDF` member.
### 3.3.2. BinnedED
This is very much one of the most important classes in the whole of OXO! This
is an N-dimensional binned event distribution. It holds a `Histogram` of the
data for the distribution, along with an `ObsSet` that defines the subset of
dimensions within the histogram that are actually observed. These objects are
typically built from MC.

With this class, you can basically do all the things you can with a `Histogram`
object.

### 3.3.3. CompositeED
Sometimes, you want to assume certain observables are independent of one another,
and still build an event distribution from them. For example, you might observe
an event's energy and radius, and instead of building a 2D `EventDistribution`
object, if you think the two observables are fairly independent then you can make
a 1Dx1D distribution instead. That's what the `CompositeED` class does: define
an event distribution as the (outer) product of some set of `EventDistribution`
objects.

This class has been defined in such a way that you can combine both `BinnedED` 
and `AnalyticED` objects, or even another `CompositeED` object!

### 3.3.4. SpectralFitDist
Unlike `AnalyticED`, `BinnedED` doesn't inherit from `FitComponent`. However,
you might want the functionality of a binned event distribution, with named
parameters for the dimensions of the underlying `Histogram`. Well, that's
what `SpectralFitDist` does!

# 4. Random number & Event generation
Generating lots of random numbers is quite useful in a variety of
places within a physics analysis. We have a single class, `Rand`,
with a bunch of static methods that act as a central place within
OXO for random number generation. There's nothing fancy here, just
a wrapper for ROOT's `TRandom3` class. You can generate random
numbers for a small number of standard distributions: uniform, Gaussian, and Poisson.

Now, with random number generation, we can generate fake datasets and event
distributions. Here are two (similar) classes that handle this:

## 4.1. DataSetGenerator
You have a whole bunch of data (e.g. MC events), and would like to generate a 
pretend sample of data, particularly useful for sensitivity studies? Well, the
`DataSetGenerator` class has you covered. Add `DataSet` objects along with their
expected rates to the generator object. The `sequential` flag for each data set
determines whether generation will use the events in each `DataSet` sequentially.
If set to false, events are selected randomly. There is also a set of `bootstrap` flags, which if set to true for a given `DataSet` will draw randomly from the
dataset with replacement: useful in a pinch if you need to have a fake dataset
of size larger than your MC!

There are two main modes of dataset generation: `ExpectedRatesDataSet()` and
`PoissonFluctuatedDataSet()`. The former will create a data set where the number
of events sampled from each `DataSet` will be equal (up to rounding) to the event
rates given. The latter will have the number of events sampled from each `DataSet`
instead as a Poisson fluctuation on the event rate: this is the most "realistic"
kind of fake data set. Note that in both modes, the output `OXSXDataSet` object
has events ordered by `DataSet` type - e.g. if a dataset was generated from
`DataSet` objects A and B, the generated dataset would have a number of A-type
events followed by a number of B-types.Both methods also have an optional parameter
which outputs by reference a vector of the number of events actually simulated
per type.

There are two other, more brute-force forms of dataset generation. The first is
`AllValidEvents()`, which just creates a `DataSet` from all of the events given,
ignoring the event rates. `AllRemainingEvents(i)`, as the name suggests, gets
back the events not used for previous event generation from the ith `DataSet`.

Finally, one can also provide cuts on what can be generated by `SetCuts()` or
`AddCut()`; the event rates the user gives then correspond to the event rates
foreach dataset *before* cuts. More on cuts in OXO [here](#61-cut).
## 4.2. BinnedEDGenerator
The `BinnedEDGenerator` works very much like `DataSetGenerator`, except the input
and output objects are no longer `DataSets` which contain sets of `Event` objects,
but instead `BinnedED` objects. After giving the `BinnedED` objects and rates for
each event type via `SetPdfs()` and `SetRates()` respectively, `ExpectedRatesED()`
and `PoissonFluctuatedED()` can generate binned event distributions as one would
expect. Becuase we are just dealing with binned histograms here, there is no need
for the additional stuff that `DataSetGenerator` has: bootstrapping, the sequential
flag, etc. are all not present here.

# 5. Data type conversions & IO
We've now introduced a bunch of OXO types that can handle data; it would be nice
to be able to both convert between them, as well as save/load them from disk.

## 5.1. DistTools
`DistTools` is a static utility class, that allows one to readily convert between
certain varieties of event distribution object.

`ToHist()` will take in a ROOT `TH1D` or `TH1D` object, or a `PDF` &
`AxisCollection` pair, and convert them to a native OXO `Histogram` object. From
this, one can then easily convert the `Histogram` into a `BinnedED` by using the
relevant `BinnedED` constructor.

For the case of 1D/2D `BinnedED` or `Histogram` objects, one can convert them
into ROOT `TH1D` or `TH2D` objects, as relevant.

## 5.2. DistFiller
Sometimes, you want to take a `DataSet` object you have, and generate
a `BinnedED` object from the events. Well, the `DistFiller::FillDist(pdf, data)` 
method does this! More precisely, the method will loop over events in the
`DataSet` object, and fill their information to the `BinnedED` given. If one
wants, a `CutCollection` and/or an `EventSystematicManager` can be given to apply
systematics and/or cuts on the events before filling. A restriction can also be
given on the number of events to fill to the `BinnedED`. Finally, one can
optionally provide a `CutLog` to store information on the cuts performed. More
information about how cuts and systematics work in OXO can be seen in Sections
[6](#61-cut) and [](), respectively.

Note: one cannot convert from a binned distribution back into a `DataSet`, as the 
event-by-event information has been lost.

## 5.3. IO
`IO` is another static utility class, this time handling saving/loading stuff
to/from files. If one wants to save a `DataSet` object to file, one can use
the `SaveDataSet()` method. Depending on the file extension given as part of
the output file path, the data will either be saved in an HDF5 file (`.h5`), or
as a `TNtuple` object within a ROOT file (`.root`). One can also directly call
the `SaveDataSetH5()` or `SaveDataSetROOT()` methods as relevant. The same is
true for `Histogram` objects via `SaveHistogram()`; note that saving a histogram
to a ROOT file is only possible if it is 1D/2D as they get saved to `TH1D`/`TH2D`
objects, as appropriate. [Note: `TH3D` objects exist! We should probably allow 
for them.]

One can similarly load `.h5` files that contain a dataset/histogram via
`LoadDataSet()`/`LoadHistogram()`. For `.root` files, a `TNtuple` object can be
loaded in directly via OXO's `ROOTNtuple` class; ROOT histogram objects must be
loaded in using the standard ROOT `TFile` approach, and then converted into a
`Histogram` via `DistTools::ToHist()` (described above).

# 6. Applying cuts to data
A major part of Particle Physics analyses is the applying of cuts an data, 
typically to maximise signal efficiency and minimise model uncertainty. Well,
OXO is able to handle applying these cuts nicely, tracking efficiencies!

## 6.1. Cut
The `Cut` class is the abstract base class for OXO that describes a single cut
on data of some variety. The defining method of this class is `PassesCut()`,
which returns a boolean based on whether a given `Event` object has passed the
cut. `Cut` objects also have a name that they can be referred to by.

OXO currently has 3 basic types of cuts, that can handle a fair number of 
situations. If you want to do anything advanced
(beyond combining cuts - we'll get to that!), then a good idea is to separately
"tag" events (e.g. an out-of-window BiPo tag), and then cut based on this tag.
Or just not do cutting with OXO - I'm not your boss.
### 6.1.1. BoolCut
`BoolCut` is arguably the simplest kind of `Cut` - `PassesCut(event)` returns
true if and only if a given observable has a specific value, e.g. one can get
whether an event has `evIndex == 0`, say. This gets defined at construction:
for this example, we could do `BoolCut b("no_retrigger_cut", "evIndex", "0");`.
### 6.1.2. LineCut
Next up,`LineCut` defines `PassesCut(event)` such that an event will pass if a
given observable is strictly greater/less than a fixed value, with the
`sidedness` flag being "upper"/"lower" determining which way.

So for example, `LineCut e_cut("energy_cut", "energy", 5.0, "lower");` defines
an energy cut to select events E > 5.0 MeV.

### 6.1.3. BoxCut
Finally, `BoxCut` has an event pass a cut if a given observable is within an
open interval, as defined by a lower & upper limit. For example, one could
make an nhit cut like so: `BoxCut n_cut("nhit_cut", "nhits_cleaned", 1000, 1500);`,
which passes if 1000 < `nhits_cleaned` < 1500.

## 6.2. CutCollection
In an analysis, one rarely applies only a single cut. What a person normally wants
to know is whether a certain event passes *all* the cuts specified: this is
what `CutCollection` is for. We add individual cuts to the collection by
`AddCut(cut)`, and then the `PassesCuts(event)` method will return true only
if all cuts are passed.

## 6.3. CutLog
A hugely useful companion to `CutCollection` is the `CutLog` class, which stores
information about the impact of cuts on data through the overloaded
`CutCollection::PassesCuts(event, log)` method. Note that before a `CutLog` object
can be successfully used in this manner, it must be constructed with the names
of the cuts being used. This can be done as easily as:`
```
CutLog log(cut_collection.GetCutNames());
```
After performing all the cutting you would like to do, there are four numbers per
cut stored, once `CalculateMeta()` is run:
  - Number of events cut after applying each cut: `GetCutCounts()`
  - Percentage of events cut after applying each cut: `GetCutPercentages()`
  - Number of events remaining after having applied all the cuts so far: `GetRemainderCounts()`
  - Percentage of events remaining after having applied all the cuts so far: `GetRemainderPercentages()`

This information can be prettily printed via `Print()`, stored as a string via
`AsString()`, or even saved to a file by `SaveAs()`.

Note: there is no in-house method I am aware of that automatically 
creates a new `DataSet` object that takes the data from some original `DataSet`
that passes all the cuts within a `CutCollection`. Here's a function that does
just that:
```
void get_data_passed_cuts(const ROOTNtuple& ntuple, const CutCollection& cuts, 
                          CutLog& cut_log, OXSXDataSet& cut_data) {
    /*
     * Iterate over events in the ntuple. If an event passes cuts,
     * add to the cut_data object, and note in cut_log.
     */
    for (size_t i = 0; i < ntuple.GetNEntries(); i++) {
        const Event event = ntuple.GetEntry(i);
        if (cuts.PassesCuts(event, cut_log)) { cut_data.AddEntry(event); }
    }
}
```

# 7. Systematics
Consideration of systematic effects is a major part of Particle Physics analyses
(and a real pain, to boot!). Fortunately, OXO allows for the handling of
systematics. There are in fact two entirely separate varieties of systematic
considered: ones applied to individual `Events` (and hence `DataSets`), versus
those applied to `BinnedED` objects.

## 7.1. EventSystematic
Let's consider systematics applied to events first: these all are derived from
the `EventSystematic` abstract base class. To be useful, one must first set the
observables upon which the systematic can affect, via `SetOutObservables()`.
Sometimes, the systematic also needs to know about the values of other
variables that aren't modified themselves: this can be set with
`SetInObservables()`.

The defining method of this class is the evaluation method, `sys(event)`. This
outputs a new `Event` object that has applied the systematic upon the input
`Event`. Note that `EventSystematic` derives from `FitComponent`, so is
able to handle any number of abstract parameters that could be needed to define
the systematic.

### 7.1.1. EventShift
`EventShift` is probably the simplest kind of systematic to apply. After
defining the observable to be shifted with `SetOutObservables()`, and setting
the shift with `SetParameter()` (technically one can also do this with 
`SetShift()`, but using `SetParameter()` allows the whole `FitComponent` to
work nicely - useful if doing fits), then one can apply the shift to an event
by evaluation. this will apply `x' = x + a`, where `x` and `x'` are the 
observable of interested before and after the transformation, and `a` is the
shift parameter.

### 7.1.2. EventScale
`EventScale` acts just like `EventShift`, except now the defining 
transformation is `x' = a*x`, using the notation defined above.

### 7.1.3. EventConvolution
TODO

### 7.1.4. EventReconvolution
TODO

## 7.2. Systematic
Unlike `EventSystematic`, the `Systematic` class is for manipulating a binned
event distribution instead of a set of events themselves. It, too is an abstract
base class which derives from `FitComponent`, so abstract paramters can be stored
within it.

Critically, we can describe a general transformation of a binned event distribution
due to some systematic effect by a linear transformation. As a result, a systematic
can be represented by a matrix (the "detector response matrix") that acts on the
`BinnedED` to produce a modified event distribution. In OXO, we use a `SparseMatrix`
class that is a wrapper for `Armadillo`'s `arma::sp_mat` class.

This response gets constructed via the `Construct()` method, which defines the 
particular kind of systematic. Any time the underlying parameters that define
the systematic get changed, this method must be called once again. The evalution
operator (`operator(binned_ed)`) can then be used to determine what impact the
detector response matrix has on the event distribution.

There are two sets of observables considered by this class: firstly, the full 
set of observables needed to describe an event distribution; and secondly the
observable(s) that the systematic actually acts upon.These can be defined by
the `SetDistributionObs()` and `SetTransformationObs()` methods, respectively.

### 7.2.1. Scale
The `Scale` systematic works like `EventScale`, except now of course the impact
of scaling a parameter on an event distribution will be approximate because the
distribution is binned.

### 7.2.2. Convolution
TODO

# Managers
In most full-blooded analyses, there are often numerous distributions, 
systematics, parameters, and fittable components, all flying around. It's 
critical to keep track of them all, especially if one is trying to perform a 
fit! That's what OXO's managers are for. Unlike most of the other classes,
these all have to do different-enough things that they don't derive from some
fundamental class; I group them here anyway as they share the same broad idea
of managing objects within OXO. These managers get used regularly within fitting.

## ParameterManager
