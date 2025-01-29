#
find_package(Catch2 3 REQUIRED)

# Define where the example macro files are
set(tests
    AxisCollectionTest.cpp
    BinAxisTest.cpp
    BinnedEDGeneratorTest.cpp
    BinnedEDManagerTest.cpp
    BinnedEDShrinkerTest.cpp
    BinnedEDTest.cpp
    BinnedNLLHTest.cpp
    ComponentManagerTest.cpp
    CompositeEDTest.cpp
    ConstraintsTest.cpp
    ConvolutionTest.cpp
    CutTest.cpp
    DataSetGeneratorTest.cpp
    DataSetIOTest.cpp
    DistToolsTest.cpp
    EDManagerTest.cpp
    EventScaleTest.cpp
    EventShiftTest.cpp
    EventSystematicManagerTest.cpp
    FitParameterTest.cpp
    GaussianFitterTest.cpp
    GaussianTest.cpp
    HistogramGetMultiDSlice.cpp
    HistogramIOTest.cpp
    HistToolsTest.cpp
    ObservableSetTest.cpp
    ParameterManagerTest.cpp
    ScaleFuncSystTest.cpp
    ScaleSystTest.cpp
    ShapeTest.cpp
    ShiftSystTest.cpp
    SpectralFitDistTest.cpp
    SquareRootScaleTest.cpp
    StatisticSumTest.cpp
    SystematicManagerTest.cpp
    VaryingCDFTest.cpp
)

add_executable(RunUnits ${tests})
target_include_directories(RunUnits PUBLIC ${OXSX_INCLUDE_DIR})
target_link_libraries(RunUnits PUBLIC oxsx)

target_link_libraries(RunUnits PRIVATE Catch2::Catch2WithMain)