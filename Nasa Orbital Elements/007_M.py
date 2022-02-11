# https://stackoverflow.com/questions/56482214/how-to-fit-an-inverse-sawtooth-function-to-a-curve-or-a-plot

import numpy, scipy, matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings
import scipy.signal


n_step_list = [-500.0, -400.0, -300.0, -200.0, -100.0, 0.0, 100.0, 200.0, 300.0, 400.0, 500.0]
value_list =  [-24.0, 73.0, 55.0, 36.0, 18.0, 0.0, -18.0, 79.0, 61.0, 43.0, 24.0]

xData = numpy.array(n_step_list)
yData = numpy.array(value_list)


# width is from scipy docs at https://www.pydoc.io/pypi/scipy-1.0.1/autoapi/signal/waveforms/index.html#signal.waveforms.sawtooth
def func(x, A, fi, offset, width):
    return A * scipy.signal.sawtooth(x / fi, width) + offset


# function for genetic algorithm to minimize (sum of squared error)
def sumOfSquaredError(parameterTuple):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    val = func(xData, *parameterTuple)
    return numpy.sum((yData - val) ** 2.0)


def generate_Initial_Parameters():
    # min and max used for bounds
    maxX = max(xData)
    minX = min(xData)
    maxY = max(yData)
    minY = min(yData)

    minData = min(minY, minX)
    maxData = max(maxY, maxX)

    parameterBounds = []
    parameterBounds.append([minData, maxData]) # search bounds for A
    parameterBounds.append([minData, maxData]) # search bounds for fi
    parameterBounds.append([minData, maxData]) # search bounds for Offset
    parameterBounds.append([0, 1]) # search bounds for width, see https://www.pydoc.io/pypi/scipy-1.0.1/autoapi/signal/waveforms/index.html#signal.waveforms.sawtooth

    # "seed" the numpy random number generator for repeatable results
    result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
    return result.x

# by default, differential_evolution completes by calling curve_fit() using parameter bounds
geneticParameters = generate_Initial_Parameters()

# now call curve_fit without passing bounds from the genetic algorithm,
# just in case the best fit parameters are aoutside those bounds
fittedParameters, pcov = curve_fit(func, xData, yData, geneticParameters)
print('Fitted parameters:', fittedParameters)
print()

modelPredictions = func(xData, *fittedParameters) 

absError = modelPredictions - yData

SE = numpy.square(absError) # squared errors
MSE = numpy.mean(SE) # mean squared errors
RMSE = numpy.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (numpy.var(absError) / numpy.var(yData))

print()
print('RMSE:', RMSE)
print('R-squared:', Rsquared)

print()


##########################################################
# graphics output section
def ModelAndScatterPlot(graphWidth, graphHeight):
    f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=100)
    axes = f.add_subplot(111)

    # first the raw data as a scatter plot
    axes.plot(xData, yData,  'D')

    # create data for the fitted equation plot
    xModel = numpy.linspace(min(xData), max(xData))
    yModel = func(xModel, *fittedParameters)

    # now the model as a line plot
    axes.plot(xModel, yModel)

    axes.set_xlabel('X Data') # X axis data label
    axes.set_ylabel('Y Data') # Y axis data label

    plt.show()
    plt.close('all') # clean up after using pyplot

graphWidth = 800
graphHeight = 600
ModelAndScatterPlot(graphWidth, graphHeight)
