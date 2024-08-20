OXSX/OXO is a signal extraction framework built for the SNO+ experiment.
It is built to be useful for general analyses within particle physics.

This document will try and explain the critical components of OXO, so you can
build your analysis code efficiently.

- [1. Objects for storing data](#1-objects-for-storing-data)
  - [1.1. Event](#11-event)
  - [1.2. ObsSet](#12-obsset)
  - [1.3. DataSet](#13-dataset)
    - [1.3.1. OXSXDataSet](#131-oxsxdataset)
    - [1.3.2. ROOTNtuple](#132-rootntuple)
    - [1.3.3. ROOTTree](#133-roottree)
    - [1.3.4. LazyOXSXDataSet](#134-lazyoxsxdataset)
- [2. Histograms, etc.](#2-histograms-etc)
  - [2.1. BinAxis](#21-binaxis)
  - [2.2. AxisCollection](#22-axiscollection)
  - [2.3. Histogram](#23-histogram)
- [3. PDFs \& Event Distributions](#3-pdfs--event-distributions)
  - [3.1. ParameterDict](#31-parameterdict)
  - [3.2. FitComponent](#32-fitcomponent)
    - [3.2.1. Function](#321-function)
      - [3.2.1.1. Heaviside](#3211-heaviside)
      - [3.2.1.2. 3.2.1.2 SquareRootScale](#3212-3212-squarerootscale)
      - [3.2.1.3. PDF](#3213-pdf)
        - [3.2.1.3.1. Gaussian](#32131-gaussian)
    - [3.2.2. ConditionalPDF](#322-conditionalpdf)
      - [3.2.2.1. JumpPDF](#3221-jumppdf)
      - [3.2.2.2. VaryingCDF](#3222-varyingcdf)
  - [3.3. EventDistribution](#33-eventdistribution)
    - [3.3.1. AnalyticED](#331-analyticed)
    - [3.3.2. BinnedED](#332-binneded)
    - [3.3.3. CompositeED](#333-compositeed)
    - [3.3.4. SpectralFitDist](#334-spectralfitdist)
- [4. Random number \& Event generation](#4-random-number--event-generation)
  - [4.1. DataSetGenerator](#41-datasetgenerator)
  - [4.2. BinnedEDGenerator](#42-binnededgenerator)
- [5. Data type conversions \& IO](#5-data-type-conversions--io)
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
    - [7.2.2. ScaleFunction](#722-scalefunction)
    - [7.2.3. 7.2.3 Shift](#723-723-shift)
    - [7.2.4. 7.2.4 Convolution](#724-724-convolution)
    - [7.2.5. 7.2.5 Shape](#725-725-shape)
- [8. Constraints](#8-constraints)
  - [8.1. QuadraticConstraint](#81-quadraticconstraint)
  - [8.2. BivariateQuadraticConstraint](#82-bivariatequadraticconstraint)
- [9. Parameter Wrangling](#9-parameter-wrangling)
  - [9.1. ComponentManager](#91-componentmanager)
  - [9.2. FitParameter](#92-fitparameter)
    - [9.2.1. DoubleParameter](#921-doubleparameter)
    - [9.2.2. ContainerParameter](#922-containerparameter)
    - [9.2.3. ParameterManager](#923-parametermanager)
    - [9.2.4. EDManager](#924-edmanager)
    - [9.2.5. BinnedEDManager](#925-binnededmanager)
- [10. Other Managers](#10-other-managers)
  - [10.1. ConstraintManager](#101-constraintmanager)
  - [10.2. EventSystematicManager](#102-eventsystematicmanager)
  - [10.3. SystematicManager](#103-systematicmanager)
  - [10.4. BinnedEDShrinker](#104-binnededshrinker)
- [11. Test Statistics](#11-test-statistics)
  - [11.1. ChiSquare](#111-chisquare)
  - [11.2. BinnedNLLH](#112-binnednllh)
  - [11.3. StatisticSum](#113-statisticsum)
- [12. Optimisers](#12-optimisers)
  - [12.1. FitResult](#121-fitresult)
  - [12.2. Optimiser](#122-optimiser)
  - [12.3. GridSearch](#123-gridsearch)
  - [12.4. Minuit](#124-minuit)
  - [12.5. Markov Chain Monte Carlo](#125-markov-chain-monte-carlo)
    - [12.5.1. MCSampler](#1251-mcsampler)
      - [12.5.1.1. MetropolisSampler](#12511-metropolissampler)
      - [12.5.1.2. HamiltonianSampler](#12512-hamiltoniansampler)
    - [12.5.2. MCMCSamples](#1252-mcmcsamples)
    - [12.5.3. MCMC](#1253-mcmc)

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
store of some observable names on their own, in a specifc order.
This class is sometimes used to define the subset of variables we actually
care about in a given event.

We could create an `ObsSet` object like:
```
const std::vector<std::string> observable_names {"energy", "radius", "ext_classifier"};
ObsSet obs(observable_names);
```

Given an `Event` `ev`, one can get the data for a particular set of 
observables
defined in an `ObsSet` `obs` via:
```
const std::vector<double> obs_data = ev.ToObsSet(obs);
```

## 1.3. DataSet

`DataSet` is just an abstract base class for some much more useful classes
we'll talk about in a moment. The point of any such object is to store a
bunch of data. The most important method here is `GetEntry(i)`, where you
get the ith event in the dataset.

Here are the derived classes of `DataSet`, which can actually be used:

### 1.3.1. OXSXDataSet

The most obvious kind of `DataSet` object in OXO: it holds a vector of `Events` 
(and the associated observable names in another vector).
In addition to `GetEntry(i)`, one can also add an event to the dataset by
`AddEntry()`. A nice feature is concatenating `OXSXDataSet` objects together
with the `+` operator overload. One can reserve memory for n events by 
`Reserve(n)`.

### 1.3.2. ROOTNtuple

A `DataSet` object derived from a ROOT `TNtuple` object that is stored in a
file. When the object is created, the `TNtuple` object is loaded. Adding
events afterwards is not possible: you're stuck with what you loaded in.

### 1.3.3. ROOTTree

Almost identical to the `ROOTNTuple` class, but loads in a ROOT `TTree` instead.
Tends to be more useful in practice than `ROOTNTuple`, because most data
one handles in particle physics gets stored in `TTree` objects.

Importantly, for this class to work the `TTree` structure must be flat, with
all branches having variables that can be cast into a double.

### 1.3.4. LazyOXSXDataSet

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
You can store them as a `ParameterDict`! etc.

This class is literally just a mapping from strings --> doubles:
```
typedef std::map<std::string, double> ParameterDict;
```

## 3.2. FitComponent
This is an abstract base class for objects that have parameters associated with
them. Anything derived from this class have some set of named parameters, each
of which is a double.

Any `FitComponent`-derived object can have parameters set & gotten (via a 
`ParameterDict` if you like!), and renamed. The object also has its own name
that can be changed.

This base class will become incredibly important when we get to test statistics
and fitting/optimisation - lots of OXO classes will be derived from `FitComponent`.

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

Frankly not really used very much.

#### 3.2.1.2. 3.2.1.2 SquareRootScale
`SquareRootScale` is a `Function` whose evaluation function is of the form
`f(x) = a*sqrt(x)`, where `a` is the "gradient" parameter. This gradient
can be set directly by `SetGradient(a)`, or indirectly using the standard
`FitComponent` interface: `SetParameter("grad", a)` etc.

This function is typcially used when modelling an energy resolution systematic: 
you expect the energy resolution to be proportional to the square root of the energy.


#### 3.2.1.3. PDF
Let's go another level deeper! a `PDF` is *another* abstract class, which is a
`Function` that also has the property of definite integration, via the 
`Integrate(mins, maxs)` method. Because a `PDF` has a well-defined integral,
it can also be also randomly sampled from: `Sample()`. See the section on 
[random number generation](#4-random-number-generation-rand) for how derived
classes implement this.

##### 3.2.1.3.1. Gaussian
For all of that effort, we can finally make some actual PDFs, in particular a
multi-dimensional Gaussian (given how useful they are!). Unsurprisingly, 
`Gaussian` derives from the `PDF` class.

This class has all of its derived features. You can get/set the
means and standard deviations of any of the dimensions by either `SetMean()` etc.,
or with the `FitComponent` `SetParameter()` method. For the latter, the name of
the mean & sigma parameters are given by `means_i` & `stddevs_i` respectively, 
where `i` here corresponds to the index of the dimension of interest.
(Note - the actual parameter names are handled not within the `Gaussian` class
itself, but within a helper class called `GaussianFitter`. Despite the name,
this latter class doesn't fit a Gaussian, but instead holds the fit parameters.)
You can evaluate the `Gaussian` at any point in the n-dimensional space, and 
even integrate within some (possibly high-dimensional) rectangular region.

As a bonus, you can obtain the cumulative density function by `Cdf()`.

### 3.2.2. ConditionalPDF
`ConditionalPDF` is an abstract class used to describe distributions 
that have conditional probabilities. These will become important later 
when dealing with convolutions: a `ConditionalPDF` object is a nice way 
of handling kernels to a convolution.

Unsurprisingly then, one of the main methods declared for this class is 
`ConditionalProbability(x,x2)`, the probability at a vector position 
`x` conditional on the position `x2`. Related to this, the `Sample(x2)` 
method should sample the PDF, conditional on the position `x2`.

This class derives from `FitComponent`, so you can have parameters associated 
with this object, as well as a name.

#### 3.2.2.1. JumpPDF
`JumpPDF` is the main derived class of `ConditionalPDF`, and is used for 
conditional PDFs of the form `f(x|x2) = f(x-x2)`. You can think of it as a 
PDF with extra translational degrees of freedom. This is typically used in situations 
where you need a PDF that must act relative to some point, e.g. during convolution.

The class holds within it an un-translated PDF `fPDF` as its base; when you calculate 
`ConditionalProbability(x,x2)` it returns `fPDF(x-x2)`. The `Sample(x2)` and 
`Integral(mins, maxes, x2)` methods behave in a similar manner.

#### 3.2.2.2. VaryingCDF
What we saw with `JumpPDF` is that you could have a PDF that, in addition to some
base parameters (which you can set/get through the `FitComponent` interface), has 
the ability to be translated. A critical feature of the `JumpPDF` class is that its 
shape is unchanged with translations.

However, you might want a translatable PDF for which its parameters have 
some proscribed dependence on the translation point: this is where `VaryingCDF` 
comes in.

Like `JumpPDF`, `VaryingCDF` derives from `ConditionalPDF`, and has many of the same 
features: there is some underlying `PDF` which can been sampled from, integrated, etc.,
with the ability to translate the PDF. The key extra feature this class posesses is that 
you can link underlying PDF parameter values to the translation position, through the 
method `SetDependance(param_name, func)`. `func` here is a `Function` object that 
must encode this dependence.

As an example, we could prepare an energy-dependent "smearer" that has a 
square-root dependence of the energy resolution with energy:
```
Gaussian gaus_energy("gaus_energy", 0, 1); // Base Gaussian kernel
SquareRootScale e_res_dependence("energy_res_scaler"); // Function that scales as a sqrt
VaryingCDF energy_smearer("energy_smearer"); // Initialise smearer object
energy_smearer.SetKernel(&gaus_energy); // Set the Gaussian as the kernel function
energy_smearer.SetDependance("stddevs_0", &e_res_dependence); // Set the Gaussian sigma param to go as e_res_dependence
```
We will see the key application of all of this when we get to convolutions.

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
One possible event distribution is a purely analytic one, given by `AnalyticED`.
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

In order to correctly set up a `BinnedED` object, you must provide to it:
 - A name
 - An `AxisCollection` that define a binning over a set of variables
 - An `ObsSet` that defines the subset of binned variables that are actually "observables"
 - (Likely) non-zero values for the contents of each bin - bin contents are zero by default.

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
If set to false, events are selected randomly. There is also a set of `bootstrap` flags, 
which if set to true for a given `DataSet` will draw randomly from the
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
events followed by a number of B-types. Both methods also have an optional parameter
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
by evaluation. This will apply `x' = x + a`, where `x` and `x'` are the 
observable of interested before and after the transformation, and `a` is the
shift parameter.

### 7.1.2. EventScale
`EventScale` acts just like `EventShift`, except now the defining 
transformation is `x' = a*x`, using the notation defined above.

### 7.1.3. EventConvolution
The `EventConvolution` class has some kernel `ConditionalPDF` object, 
and "smears" the observable of interest by sampling from the kernel, using 
the `ConditionalPDF::Sample()` method. As an example, you could smear the 
energies of events with a `Gaussian`.

### 7.1.4. EventReconvolution
`EventReconvolution` behaves very differently to `EventConvolution`. For this 
latter class, it is assumed that you already have variables in your event 
corresponding to both pre- and post-smear. The truth and reconstructed energy 
of an event are classic examples. Instead of doing more convolutions, this class 
just modifies the post-smear value linearly relative to the pre-smear value: 
`post_smear_new = pre_smear + fCorrection*(post_smear_old - pre_smear)`. 

This class has one parameter, `correction`, which defines the smearing scale. 
If `correction = 1`, nothing changes. If `correction = 0`, all smearing is removed 
and the post-smear values will equal that of pre-smear.

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
observable(s) that the systematic actually acts upon. These can be defined by
the `SetDistributionObs()` and `SetTransformationObs()` methods, respectively.

### 7.2.1. Scale
The `Scale` systematic works like `EventScale`, except now of course the impact
of scaling a parameter on an event distribution will be approximate because the
distribution is binned.

### 7.2.2. ScaleFunction
The `ScaleFunction` class behaves like `Scale`, but allows for the more complex  
scenario where you want the scaling to be dependent in some way on the observable 
values. As an example, you might want to have a non-linear energy scale systematic.

In order to do this, the user must provide to the `ScaleFunction` object a 
`std::function<double(const ParameterDict&,const double&)>` object which, 
given the value to be scaled and the `ParameterDict` of `ScaleFunction` 
parameters, returns a scaled value. This kind of function has been typedefed as
`ScaleFunc`.

### 7.2.3. 7.2.3 Shift
The `Shift` systematic works like `EventShift`, shifting a distribution by some 
constant amount.

### 7.2.4. 7.2.4 Convolution
The `Convolution` systematic works similar to `EventConvolution`, but now smears 
each bin's contents via calls to a `ConditionalPDF::Integral()` method. Requires 
a `PDF` or `ConditionalPDF` to be given for the smearing kernel.

### 7.2.5. 7.2.5 Shape
The `Shape` systematic applies a bin-by-bin scaling of its contents, effectively 
modifying the overall shape of the distribution. The shape function must be 
provided by the user as an `std::function<double(const ParameterDict&,const double&)>` 
object, typedefed to `ShapeFunction`.


# 8. Constraints
A standard component of test statistics are contraint penalty terms. In OXO, there 
is currently no fully-general way of handling constraints. There are currently two 
kinds of constraint classes that handle many standard situations.

## 8.1. QuadraticConstraint
The `QuadraticConstraint` class handles the very standard situation of wanting a 
penalty function on a single parameter of the form:
```
f(x) = (x - mu)^2/(2*sigma^2)
```
or:
```
f(x) = (x - mu)^2/(2*s_low^2)  if x <  mu,
     = (x - mu)^2/(2*s_high^2) if x >= mu.
```
The first case is a standard quadratic penalty with a mean and width parameter; 
the second is an asymmetric quadratic penalty with different widths depending 
on whether the function is evaluated above or below the mean param value.

The constraint parameters are set up with the `SetConstraint()` methods; the 
penalty is evaluated for a given value with `Evaluate(x)`.

## 8.2. BivariateQuadraticConstraint
The `QuadraticConstraint` class is limited in part because of not allowing for 
penalty terms to account for correlations between multiple parameters. Currently 
in OXO, we only have one class that handles a correlated constraint: `BivariateQuadraticConstraint`. 
This class allows one to evaluate a constraint of two parameters, where there is 
some amount of correlation between them (correlation must be < 1 though).

As expected, the penalty term can be evaluated with the `Evaluate(x1, x2)` method.

# 9. Parameter Wrangling
We've seen that a large fraction of classes within OXO inherit from `FitComponent`, 
and so they have parameters that can be varied. We're going to see when it comes to 
test statistics and optimisations that knowing what the status of all the parameters 
are for all the different bits of the analysis at the same time, and having one 
interface for changing them, will be really valuable.

There are in fact two separate ways of doing this within OXO, both which get used. 
One approach is relatively straightforward; the other is somewhat more complex.

Note that the classes in this section are usually back-end classes: it will be rare 
that you should have to interact with these classes directly.

## 9.1. ComponentManager
The `ComponentManager` class is the simple, and somewhat logical conclusion to having 
loads of objects that derive from `FitComponent`. It stores a vector of pointers to 
`FitComponent` object that you want to keep track of. To add an object to be tracked, 
you can simply use the `AddComponent(obj)` method.

With this class, you can then get and set all of the parameters in all of the tracked 
`FitComponent` objects simultaneously, via `GetParameters()` and `SetParameters(dict)`. 
There are a few other functions that allow you to act over the whole collection of 
tracked objects simultaneously. Neat!

## 9.2. FitParameter
The first method works fine, as long as you're dealing with `FitComponent` objects. 
Some OXO objects don't inherit from `FitComponent` though, e.g. `BinnedED`, so this won't 
work in general!

To solve this problem, there is a second approach to handle parameters. It is more 
complicated, but in the end we will see that it grants `FitComponent` status to 
collections of `BinnedED` objects, for example.

To start with, the `FitParameter` abstract base class defines the concept of a 
single parameter, which merely has the `Get()` and `Set()` methods.

NOTE: `FitComponent` and `FitParameter` are totally different things! Both are abstract 
base classes, but `FitParameter` describes a single parameter, whereas `FitComponent` is 
any object that has fittable components. Sorry.

### 9.2.1. DoubleParameter
The simplest implementation of `FitParameter`, `DoubleParameter` is just a class 
that stores a reference to a single double, which can be `Get()` and `Set(x)`.

### 9.2.2. ContainerParameter
The `ContainerParameter` class allows you to create a `FitParameter` object for 
any double that's within any reasonably-defined container object, such as a `std::vector`. 
All you have to do is provide the container and the index of the parameter within that 
container:
```
std::vector<double> vals = {2.5, 6.0, 3.2};
// Track the second value in the above object:
ContainerParameter cp(vals, 1);
double x = cp.Get(); // x = 6.0;
cp.Set(3.0); // now vals = {2.5, 3.0, 3.2}
```

### 9.2.3. ParameterManager
Using `DoubleParameter` and `ContainerParameter`, we can now endow basically any double 
variable in our code as a `FitParameter` object. In order to corall all of these together, 
and generally provide some quality-of-life features, the `ParameterManager` is used.

The `ParameterManager` class merely contains a mapping of parameter names, and their 
associated `FitParameter` object pointers. Much like `ComponentManager`, you can get 
and set any parameters that are being tracked by the manager object, using the usual 
`GetParameters()`, `SetParameters()` methods etc.

In order to track parameters with this class, you could create your own `DoubleParameter` 
or `ContainerParameter` object, as appropriate, and then run the `ParameterManager::Add(fit_param, param_name)` 
method for each parameter. This is likely to get tedious! Instead, you can get the class 
to do all of the `FitParameter` object creation for you using the `AddDouble()` and `AddContainer()` 
methods.

### 9.2.4. EDManager
So far, all of the `FitParameter` stuff has been pretty hypothetical: now to provide a couple 
of concrete uses for this second approach, as well as a way to bridge from `FitParameter` 
world to `FitComponent` world. After all, we eventually want to have /all/ of the parameters 
we care about handleable through a single `ComponentManager` object.

The `EDManager` class is an example of such a bridge. It contains a collection of 
pointers to `EventDistribution` objects, as well as a vector with their associated 
normalisations. A `ParameterManager` object within the class then keeps track of all 
the normalisation values.

The key thing that the `EDManager` class can do is the `Probability(event)` method, 
which returns the sum of all probabilities of a given `Event` for each 
`EventDistribution` object being tracked, weighted by the normalisations.

The `EDManager` class inherits from `FitComponent`, with the `FitComponent` 
interface simply referring to that of the `ParameterManager`. By consequence, 
it can be added to a central `ComponentManager` object.

### 9.2.5. BinnedEDManager
The `BinnedEDManager` class works a bit like `EDManager`, but handles the more specific 
case of a collection of `BinnedED` objects. It also contains a `ParameterManager` to 
track normalisations; it inherits from `FitComponent` similarly.

This class also allows for some advanced things to be done with the binned distributions. 
A collection of systematics can be applied to the `BinnedED` objects through 
the `ApplySystematics(sys_man)` method, where `sys_man` is an object of 
the `SystematicManager` class (more on this class in a bit). The `fOriginalPdfs` 
member contains all of the PDFs before systematics have been applied; 
`fWorkingPdfs` after. The `Probability()` and `BinProbability()` methods 
always refer to the latter collection.

The `ApplyShrink(shrinker)` method allows one to "shrink" the `BinnedED` 
objects via a `BinnedEDShrinker` object (more of this class in a bit). Once 
again, this impacts `fWorkingPdfs`.

An annoying consequence of applying some systematics is that the normalisations 
of `BinnedED`s can change! For example, applying a shift to a distribution could 
shift some of the events out beyond the range of the distribution. This becomes 
really important with test statistics, where we often need to know what the 
PDF normalisations actually are! We have to then keep separate in these situations 
the original normalisations before systematics have been applied, `fFittableNorms`, 
and the post-systematic normalisations, `fNormalisations`.

The `NormFittingStatus` enum, with choices `FALSE`, `DIRECT`, and `INDIRECT`, defines 
how the normalisation of each `BinnedED` object gets handled. The default mode, 
`DIRECT`, behaves in the naive way: the PDF's normalisation is a tracked `FitParameter`, 
and if systematics/PDF shrinking should impact the normalisation, the change is ignored.

For most situations, the `INDIRECT` mode is more appropriate: the pre-systematic 
normalisation is a `FitParameter`, and there is a separate pos-systematic/shrinking 
normalisation that can be acquired through `GetNormalisations()`.

Finally, you might have a PDF which you want to allow a systematic object to control the 
normalisation directly, and have the pre-systematic normalisation variable within 
the `BinnedEDManager` turned off. This is the `FALSE` mode.


# 10. Other Managers
In the previous section, a bunch of "Manager" classes were explained. Well, it turns out 
there are even more of them in OXO!

## 10.1. ConstraintManager
The `ConstraintManager` class stores a collection of constraint objects - either 
`QuadraticConstraint` or `BivariateQuadraticConstraint` at the moment - and allows the 
user to centrally control the constraint parameters via `SetConstraint()` methods. 
It also allows users to simultaneously `Evaluate(params)` them in one call, and it will 
sum the totals of all `Evaluate()` calls for each tracked constraint.

## 10.2. EventSystematicManager
If you have a bunch of `EventSystematic` objects, and want to deal with them in one place, 
then the `EventSystematicManager` is the class for you. It stores a collection of pointers 
to `EventSystematic` objects, which can be added to with the `Add(ev_sys)` method.

The key method here is `ApplySystematics(event)`: an `Event` object will have each 
systematic that is being tracked applied to it. Neat!

## 10.3. SystematicManager
Similar to `EventSystematicManager`, the `SystematicManager` class stores a collection 
of `Systematic` object pointers. It also allows for `Systematic` objects to be tracked 
through `Add(sys)`. There are some key advancements in what this class can do, however.

Firstly, given fixed systematic parameters, a `Systematic` object will always generate 
the same response matrix via running `Systematic::Contruct()`, followed by 
`Systematic::GetResponse()`. Therefore, we can cache the response matrix for each 
systematic that is being tracked all at once when we run `SystematicManager::Construct()`. 
Even better, at the same time we can pre-multiply all of the response matrices together 
to create a "total response" matrix. It is this total response matrix which can be applied 
to `BinnedED` objects using the `DistortEDs(orig_EDs, smeared_EDs)` method.

As a bonus, the post-systematics normalisations for each of the `BinnedED` objects can be 
obtained by providing an extra `std::vector<double>* norms` object: the normalisations are 
provided by reference.

The default behaviour of this class is to apply all systematics to any `BinnedED` objects in 
an identical manner. This might not be the behaviour you want, though! It's very normal 
to have an analysis where certain systematics only apply to specific "groups" of PDFs. 
Fortunately, the `SystematicManager` class also allows for this functionality. 
You can add a `BinnedED` object to a group (or multiple groups!) through the 
`AddDist(pdf, group_name)` methods. You can then associate a systematic with a group when 
adding it to be tracked by the manager: `Add(sys, group_name)`. If a systematic is added 
without a group, then it is "global", and is applied to all PDFs.

## 10.4. BinnedEDShrinker
The `BinnedEDShrinker` isn't really a "manager", but it's worth talking about at this point 
anyways. When applying systematics to PDFs, problems can rapidly arise near the edges of 
the PDFs. In order to avoid this issue, it is generally wise to put "buffer" bins around the 
edges of the PDFs. The buffer bins continue to store PDF information as usual, so that 
systematics can still be applied in the region, but are then sacrificed during the final 
test statistic calculation.

This "sacrificing" off certain bins of the edge of a `BinnedED` object is what the 
`BinnedEDShrinker` class is for. The `SetBuffer(axis_name, num_low, num_high)` sets 
the number of bins within a given axis at the bottom and top to declare as buffers. 
If you have multiple axes that have systematics you want buffers for, that is okay. 
Once the buffers have been declared, a mapping from bins before and after the 
shrinking has been applied needs to be cached: this is done with `SetBinMap(dist)`. 
Then, the sacrificing can be performed speedily with the `ShrinkDist(dist)` method.

The default behaviour is for the contents of the bin buffers to be lost during shrinking. 
However, you may want those contents to be put into the underflow/overflow bins: if so, 
use the `SetUsingOverflows(true)` method.


# 11. Test Statistics
A key part of most analyses is the calculation of some kind of test statistic, be it 
a chi-squared, log-likelihood, or otherwise. We'll see later how optimisers are 
based around repeatedly evaluting a test statistic at different points in the 
parameter space.

In OXO, `TestStatistic` is an abstract base class for test statistics. It requires 
any test statistic to have certain methods for setting/getting fit parameters, 
such as `SetParameters(params)` and `GetParameters()`, a `RegisterFitComponents()` 
method to set up all the fit components within the relevant `ComponentManager`, 
and most importantly an `Evaluate()` method to get the test statistic value for 
a given set of parameter values. The `RegisterFitComponents()` method must be 
run once before any evaluations are performed.

## 11.1. ChiSquare
The `ChiSqaure` class is a `TestStatistic` that calculates the chi-squared value 
for a collection of `BinnedED` objects when compared to a data distribution. This 
chi-square value corresponds to `sum((obs_i - mc_i)^2/mc_i)`, where `obs_i` is the ith 
data bin contents, and `mc_i` is the summed (and normalisation-weighted) PDF bin contents.

WARNING: this class is currently broken! It requires `BinnedED` PDFs to be provided to it 
for comparison to data, but has no way of adding a PDF to the class.

## 11.2. BinnedNLLH
The binned negative log-likelihood is a test statistic used ubiquitously within 
experimental particle physics, and is implemented in OXO with the `BinnedNLLH` 
class. This class has many features, allowing fairly complicated analyses to 
be performed.

This class contains many of the classes already discussed in the rest of this document. 
Let's go through the functionality.

The first thing that needs to be provided to this class is some data! If you already have 
a pre-binned data distribution in `BinnedED` form, it can be added via `SetDataDist(data_dist)`. 
Alternatively, you may want to provide an actual `DataSet` object via `SetDataSet(data)`, 
and the class will do the binning of the data automatically. If you do the latter 
approach, you can also get it to apply a series of cuts via `AddCut(cut)` or `SetCuts(cuts)`,  
where `cut` and `cuts` are `Cut` and `CutColection` objects, respectively. The class contains 
a `CutLog` object to see the results of these cuts being applied to the data.

The second required thing for `BinnedNLLH` to work are any number of `BinnedED` objects. 
These can be added with the `AddPdf()` and `AddPdfs()` methods; there are numerous overloaded 
versions. All of the `BinnedED`s are stored within a `BinnedEDManager` object; as a result 
one of the optional arguments here is to specify what the `NormFittingStatus` values are for 
each PDF (see the subsection on `BinnedEDManager` for more details). If you really want to, 
you can create your own `BinnedEDManager` object and then give it to `BinnedNLLH` via the 
`SetPdfManager(pdf_manager)` method.

Once you've provided this information, you can run the `RegisterFitComponents()` so that the
`ComponentManager` object can know about all of the relevant PDF normalisation fit parameters. 
By doing this, `BinnedNLLH` can now use the rest of it `TestStatistic` interface, which will 
be necessary when it comes to using `Optimisers` later on. If not using an optimiser, you'll 
want to set the fit parameters manually via `SetParameters(params)`.

This is the minimum you need to do for the most important method, `Evaluate()` to 
actually work. This will calculate the extended binned nagative log-likelihood, between 
the data and normalisation-weighted sum of the PDFs. Nice!

There are a whole bunch more advanced features present, though. Constraints on parameters 
in the fit can be added via one of the `SetConstraint()` methods; this adds either a 
`QuadraticConstraint` or `BivariateQuadraticConstraint` object to the `BinnedNLLH`'s 
`ConstraintManager` object.

The other key feature is support for systematics. Systematics can be added through 
the `AddSystematic(sys)` and `AddSystematics(syses)` methods, along with their overloads. 
These systematics are stored within a `SystematicManager` object, and so there is support 
for groups of PDFs for which different systematics can be applied: these can be defined in 
the `AddSystematic()` and `AddPdf()` methods appropriately. If you really want to, you 
can manually create your own `SystematicManager` object and give it to `BinnedNLLH` via 
the `SetSystematicManager(sys_man)` method.

In order for the systematics 
to be correctly applied to the PDFs during an `Evaluate()` call, `RegisterFitComponents()` 
must first be run to register the systematic fit parameters within the `ComponentManager`, 
and the systematic fit parameters set to some values via `SetParameters(params)`.

If using systematics here, you'll want to have thought carefully about what the 
appropriate `NormFittingStatus` values are for each of your PDFs. In addition, 
you should decide if adding buffer bins to the edges of certain axes is appropriate. 
The class contains a `BinnedEDShrinker` object that allows you to handle this, via 
the `SetBuffer(ax, n_lo, n_hi)` and `SetBufferAsOverflow(bool)` methods.

A source of uncertainty in an analysis that is rarely-discussed but often worth 
thinking about is the number of events simulated to build a given PDF. If the 
MC statistics was low, then that should increase the uncertainty and hence 
broaden the log-likelihood contour. One method for handling this has been 
developed by Barlow & Beeston in https://arxiv.org/abs/1103.0354. This 
has been implemented into the `BinnedNLLH` class, and can be turned on with 
`SetBarlowBeeston(true)` (it is turned off by default). Using this approach 
will require the MC statistics for each PDF to be given when using the 
`AddPdf()` or `AddPdfs()` methods.

Finally, because all of the above features can make the evaluation calculation a 
bit complicated at times, there is an option to turn on a "debug mode" via 
`SetDebugMode(true)`. This will add copious print statements when running 
`Evaluate()`, given full details about how the calculation is performed. 
This can be very useful when you're confused about what has gone wrong in your fit!

## 11.3. StatisticSum
What if you're in a situation where you have multiple different datasets that you 
want to fit simultaneously with one combined test statistic? This is where the 
`StatisticSum` class comes in. For each dataset, set up the relevant `TestStatistic` 
object, and then add each of them to `StatisticSum` via `AddStat(test_stat)`. You 
can even add them with the "+" operator overloads, e.g.:
```
ChiSquare t1;
BinnedNLLH t2;
StatisticSum t_comb = t1 + t2; // Lovely!
```

Once all the test statistics have been added, you can use the single `TestStatistic` 
interface within the `StatisticSum` object to handle things: `RegisterFitComponents()`, 
`SetParameters(params)`, etc. The `Evaluate()` method in particular will return 
the sum of all the registered `TestStatistic`'s `Evaluate()` methods.


# 12. Optimisers
Test statistics are only so useful on their own. At some point, we want to map out the 
fit parameter space, or find the point in parameter space that optimises for the
maximum/minimum test statistic value. OXO currently has features for three different 
kinds of optimiser: basic grid search, Minuit, and Markov Chain Monte Carlo. Each 
behave very differently.

## 12.1. FitResult
Before we go anywhere, we need a way of providing the final result of a fit in some way: 
this is what the `FitResult` class does. This is a glorified struct, storing the 
best-fit result as a `ParameterDict` that can be obtained through `GetBestFit()`. 
The object may have additional information, depending on the optimiser used: the 
covariance matrix, the value of the test statistic at the best-fit value, whether the 
fit result is valid, and even `Histogram` objects of the test statistic space.

The only non-getter/setter methods for this class are `AsString()`, which provides a 
nicely-formatted string with the best fit results, and `SaveAs(filename)`, which 
writes that string to a file.

## 12.2. Optimiser
`Optimiser` is an abstract base class that defines what an optimiser is in OXO. 
The only requirement for any `Optimiser` is the existence of an `Optimise(test_stat)` 
method, which takes a `TestStatistic` pointer and returns a `FitResult` object.

## 12.3. GridSearch
Arguably the simplest kind of optimisation algorithm is that of a grid search: just 
calculate the test statistic value at every point in a regularly-spaced n-dimensional 
grid of positions in the n-dimensional parameter space, and see which one gives you 
the maximal value. This is encoded within the `GridSearch` class.

Beyond the setup of the `TestStatistic` object that `GridSearch` needs to optimise over, 
the user must specify the minima, maxima, and step sizes of the search grid along all 
off the fit parameters, via `SetMinima(minima)`, `SetMaxima(maxima)`, and 
`SetStepSizes(steps)`.

By default, `GridSearch` will find the minimum value of the test statistic. However, 
you can flip this to find the maxmimum via `SetMaximising(true)`.

The only information provided by `GridSearch` in the `FitResult` object is the 
best-fit point, and the associated test statistic value.

## 12.4. Minuit
The Minuit optimiser is a classic set of optimsiation algorithms used widely in particle 
physics. If you've set the optimisation hyper-parameters up sufficiently well, 
then Minuit is often great at optimising over the parameter space in a robust way. 
The OXO `Minuit` class is just a wrapper for the C++ implementation of these 
algorithms that are stored within the `ROOT` external library provided by the user, 
so that the interface works within the OXO framework. Note that we do not try 
to provide full documentation about the underlying Minuit algroithms, just the 
OXO interface. (There is plenty of documentation for this on the internet!)

In order to set up `Minuit` correctly, you must first choose which minimisation 
algorithm to use: "Migrad" (the default), "Minimize", or "Simplex". This can be 
set with `SetMethod(method_str)`. Then, you can provide initial values, initial 
errors, minima, and maxima for the parameters through the relevant setter methods. 
Along with this, the maximum number of calls to the test statistic evaluation method, 
and the optimiser tolerance, should also be provided.

The `FitResult` generated with this optimiser provides the best-fit point and its 
test statistic value, whether the fit was valid, and the covariance matrix. In order 
to make sure the covariance matrix is correct, the difference in test statistic value 
that defines a 1 sigma error must be given, via `SetUpperContourEdge(num)`. If using 
`ChiSquare` this is 1, and if using `BinnedNLLH` it is 0.5.

Parameters can be fixed in the fit via the `Fix(param_name)` method. The fitter can 
also be set to maximise instead of minimise via `SetMaximising(true)`.

Under the hood of the `Minuit` class, the `TestStatistic` object is stored within a 
`MinuitFCN` class, which inherits from Minuit's `FCNBase` class. Minuit needs such 
a class to make the evaluation calls on. Fortunately, no direct interaction with 
this class should ever be needed.

## 12.5. Markov Chain Monte Carlo
The final `Optimiser` currently implemented within OXO is that of Markov Chain 
Monte Carlo (MCMC). Unlike other optimisation algorithms, 
MCMC doesn't try and find the position in the parameter space with the maximum/minimum 
value of the test statistic. Instead, it takes random steps through the parameter 
space, such that the overall sampled distribution of parameters tends towards 
the posterior density distribution of the parameter space, assuming the test statistic 
corresponds to the sum of the log-likelihood and prior distributions of the space.

### 12.5.1. MCSampler
To get MCMC to work, a critical part is the choice of algorithm for making random 
samples of the parameter space. `MCSampler` is an abstract base class for these 
sampling algorithms, of which two have been implemented in OXO: the "Metropolis" 
algorithm, and "Hamiltonian" Monte Carlo.

Any `MCSampler` class must provide a `Draw(current_step)` method, which proposes 
a new position in the parameter space to try and step to, given the current step. 
There must also be the method `CorrectAccParam()`, which provides a correction 
factor to the eventual acceptance probability for the proposed step.

#### 12.5.1.1. MetropolisSampler
The "Metropolis" algorithm is implemented within the `MetropolisSampler` class, 
which derives from `MCSampler`. In many ways, it is one of the simplest possible 
MCMC algorithms. The `Draw(current_params)` method here consists of displacing 
the current position by sampling according to a Gaussian for each parameter. 
The Gaussian is centered around the current position, with the widths of the 
Gaussian in each parameter set beforehand via `SetSigmas(sigmas)`.

No correction to the eventual acceptance probability is needed with this algorithm, 
so `CorrectAccParam()` always returns zero (the log of the correction to the 
acceptance probability).

Note that this is only a simplified version of the somewhat more famous 
"Metropolis-Hastings" algorithm. In the latter, you can deliberately have the 
proposal distribution not be symmetric with respect to the current and proposed 
step, which can be really useful e.g. at parameter boudaries. We do not 
have this implemented in OXO currently; you can only use the symmetric Gaussian 
proposal distribution. This can mean that for an analysis where the MCMC chain 
is often trying to be near a high-dimensional "corner" of the parameter space, 
the acceptance probability could be quite low.

#### 12.5.1.2. HamiltonianSampler
A more complex MCMC sampling algorithm is that of "Hamiltonian Monte Carlo" (HMC), 
defined in OXO with the `HamiltonianSampler`. Here, 
the test statistic is thought of defining a "potential" in the parameter space, 
and at every step in the chain a mini metaphorical physics "simulation" is performed, 
of a mass thrown in some random direction, moving through the potential field. 
The end position of this simulation defines the proposed point in parameter space 
returned by the `Draw(current_params)` method.

Unlike the Metropolis algorithm, HMC requires knowledge of the gradient of the 
test statistic value in order to perform the "simulation". This is 
estimated numerically with a custom `Gradient` class; it also means that 
`HamiltonianSampler` must be given the `TestStatistic` object on construction.

This class must be provided with "mass" hyper-parameters for each fit parameter 
via the `SetMasses(masses)` method. The square roots of these masses 
determine the widths of Gaussians that are sampled for the initial momentum 
components along the axis of each fit parameter. On construction of the 
sampler, the number of steps to numerically simulate, as well as a global 
relative step size parameter "epsilon" for the simulation must also be given.

In order for this sampler to achieve detailed balance (critical for MCMC), the 
(log of the) acceptance probability `CorrectAccParam()` must be corrected by the 
overall change in the simulated particle's kinetic energy.

The numerical simulation of this "particle" is handled within the `LeapFrog` 
static class methods. Importantly, it includes handling of fit parameter spaces 
with boundaries: the "particles" reflect off of the boundary "walls". In 
doing so, this elevates the algorithm to "Reflective HMC" (RHMC). RHMC is 
known to be far more efficient at proposing parameter positions than regular HMC. 
For this to work, the `HamiltonianSampler` class must be provided the fit 
parameters' minimum and maxmimum values, via `SetMinima(minima)` and `SetMaxima(maxima)`.

Note: don't get confused between the steps that the HMC algorithm takes, and the 
overall MCMC alogrithm!

### 12.5.2. MCMCSamples
At the other end of the MCMC process, we would like a away to store information about 
what has happened when we ran the chain. Most of the statistical information is stored 
in the overall distribution of the steps in 
the chain. The `MCMCSamples` class stores this information, along with sundry metadata 
about the MCMC fit.

Most of the methods in `MCMCSamples` are used internally by the `MCMC` class. The most 
important class that can be used once the `MCMCSamples` class has already been filled 
is `GetChain()`, which will return a ROOT `TTree` with information about the MCMC 
chain.

### 12.5.3. MCMC
Finally we reach the actual `MCMC` class! This is the actual `Optimiser` class, with 
the appropriate `Optimise(test_stat)` method to run the MCMC. You will need to set 
some things up first, though.

The first requirement is to provide the `MCMC` object with a `MCSampler` object at 
creation, defining the choice of MCMC algorithm to be used (Metroplis or RHMC 
currently). Then, the number of steps in the chain to run, the initial step 
position in the fit parameter space, and the minimum and maximum value of all the 
fit parameters must all be provided through the relevant setter methods.

There are a bunch of additional settings for the MCMC that can be modified. The 
sign of the test statistic can be flipped with `SetFlipSign(true)`: this is needed 
if using the `BinnedNLLH` test statistic. Also needed in that case is the 
`SetTestStatLogged(true)` method, which says that the test statistic corresponds to 
the logarithm of the probability density. 

An `MCMCSamples` object can be returned after the `Optimise()` method has been run. 
If you would like that object to contain a `TTree` with the step positions of the 
chain, make sure to do `MCMC::SetSaveChain(true)` before fitting. Similarly, an 
n-dimensional `Histogram` object that stores the binned posterior density over 
the whole fit parameter space can be made with `SetSaveFullHistogram(true)`; you'll 
have to provide the binning for the axes via `SetHistogramAxes(axes)`. If set to 
false, 1D and 2D `Histograms` will still be made for all of the parameter and 
parameter pair slices. If you would 
like to not save the first N steps of the chain to the `MCMCSamples` object, you can 
use the `SetBurnIn(N)` method. If you would like only every kth event to be saved, 
then you can use `SetThinFactor(k)`.

During the running of the MCMC, the auto-correlation of the test statistic is tracked 
with the `AutoCorrelationCalc` class that is stored within the `MCMCSamples` object. 

Finally, let's finish with an example to show how this comes together. Suppose you 
have already correctly set up a `BinnedNLLH` object, `lh_function`, and you 
would like to run an MCMC over it using the RHMC algorithm. Here's the code:
```
// First - set up sampling algorithm
const double epsilon = 0.1;
const int n_steps_hmc = 100;
HamiltonianSampler hm_sampler(lh_function, epsilon, n_steps_hmc);
hm_sampler.SetMasses(masses); // Requires the masses ParameterDict to be set up!
hm_sampler.SetMinima(minima); // Requires the minima ParameterDict to be set up!
hm_sampler.SetMaxima(maxima); // Requires the maxima ParameterDict to be set up!

// Now, set up the MCMC object itself
MCMC mcmc(hm_sampler);
const unsigned n_steps = 1000000;
mcmc.SetMaxIter(n_steps);

mcmc.SetFlipSign(true);
mcmc.SetTestStatLogged(true);

mcmc.SetMinima(minima);
mcmc.SetMaxima(maxima);
mcmc.SetInitialTrial(init_vals); // Requires the init_vals ParameterDict to be set up!

mcmc.SetHistogramAxes(param_axes); // Requires the param_axes AxisCollection to be set up!

mcmc.SetSaveChain(true);
const unsigned burn_in = 10000;
mcmc.SetBurnIn(brun_in);

// Okay, we're good to go! Run the MCMC!
const FitResult& best_fit = mcmc.Optimise(lh_function);
// Get the detailed chain information
const MCMCSamples mcmc_samples = mcmc.GetSamples();
TTree* chain_tree = mcmc_samples.GetChain();
const std::map<std::string, Histogram>& all_1d_hists = mcmc_samples.Get1DProjections();
```