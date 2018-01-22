from __future__ import division,print_function
from subprocess import call
import pandas as pd
import numpy as np
import scipy.io
from os import remove,getcwd,listdir
from os.path import join,dirname,isfile


def octave_file():

    '''Check whether Octave program simulation.m exists and,if it doesn't,create it.'''
    if not isfile('simulation.m'):

        lines = ['% simulation.m\n',
                    'function [] = simulation;\n',
                    '\n',
                    'format long;\n',
                    '\n',
                    '% Assign values from arguments\n',
                    'arg_list = argv();\n',
                    'filename = arg_list{1};\n',
                    'moduleDir = arg_list{2};\n',
                    '\n',
                    'addpath(moduleDir)\n',
                    '\n',
                    '% Load shocks\n',
                    'load shocks.mat;\n',
                    '\n',
                    '% Simulate model\n',
                    'data = {};\n',
                    'for i =1:numel(shocks);\n',
                    '\tdata{i} = dynare_simul(filename,shocks{i});\n',
                    '\n',
                    '% Export simulated data\n',
                    "save('-mat','simulated_data.mat','data');\n",
                    '\n',
                    'end;']

        with open('simulation.m','w') as newfile:
            newfile.writelines(lines)

octave_file()

class model:

    '''This program defines a class -- dynarehelper.model -- with associated methods for solving and 
    simulating dynamic stochastic general equilibrium (DSGE) models using Dynare++ and Octave.'''

    def __init__(self,parameters=None,shockNames=None,varNames=None,equations=None,steadyStates=None,covMat=None,order=1,filename='model.mod'):
        
        '''Initializing an instance dynarehelper.model requires values for the following varNames:

        Args:
            parameters:     (Pandas Series) contains the parameters of the model.
            shockNames:     (list) A list of strings with the names of the exogenous shocks
            varNames:       (list) A list of strings with the names of the endogenous variables.
            equations:      (list) A list with string elements that correspond to equation lines 
                                for a Dynare++ file. No semicolons.
            steadyStates:   (list,Numpy array,Pandas Series) Initial values for the Dynare++
                                nonlinear solver. Can be a Pandas Series with index equal to varNames,
                                a numpy array with values ordered to conform with varNames,or a list
                                strings with equations to be evaluated by the solver
            covMat:         (Numpy array) Shock covariance matrix
            order:          (int) Order of approximation
            filename:       (str) nome of .mod file to be created for Dynare++ computation. If string
                                doesn't end with .mod,then .mod is appended

        Returns:
            None

        Attributes:
            parameters:     (Pandas Series) 
            shockNames:     (list)
            varNames:       (list) 
            equations:      (list)
            steadyStates:   (list,Numpy array,Pandas Series)
            covMat:         (Numpy array)
            order:          (int)
            filename:       (str)
            moduleDir:      (str)

        '''

        # Assign attribute values
        self.parameters = parameters
        self.shockNames = shockNames
        self.varNames = varNames
        self.equations = equations
        self.steadyStates = steadyStates
        self.covMat = covMat
        self.order = order
        if filename.endswith('.mod'):
            self.filename = filename
        else:
            self.filename = filename+'.mod'

        self.filename = filename
        self.moduleDir = dirname(__file__)

        # Create Dynare++ .mod file
        create_dynare_file(parameters,shockNames,varNames,equations,steadyStates,covMat,order,filename)


    def approximate(self,per=0,dropFirst=0,sim=0,centralize=True,order= None,otherOptions=None):

        '''Use Dynare++ to approximate the decision rule for the nonlinear model

        Args:
            per:            (int) Number of periods to simulate. Default: 0
            dropFirst:      (int) Number of initial simulated periods to discard. Default: 0
            sim:            (int) Number of simulation periods. Default: 0
            centralize:     (bool)=False. Whether to subtract the constant term of the decision rule (not 
                                necessarily steady state for higher order approximations). Default: True
            order:          (int) Order of approximation. Overrides self.order. Default: None
            otherOptions:   (list) additional arguments to be passed to Dynare++.  Default: None

        Returns:
            None

        Attributes:
            ss:         (Pandas Series) Steady states computed by Dynare++ nonlinear solver
            shell_cmd:  (str) shell command used to run Dynare++
            stateNames: (list) list of state variables
            shockNames: (list) list of shocks
            locations:  (dict) identifies the ordering (starting from 0) of the variables in the model
                            as determined by Dynare++

        '''

        # Approximate the decision rule
        results = run_dynare(self.filename,per=per,dropFirst=None,sim=sim,centralize=centralize,order=order,otherOptions=otherOptions)

        # Assign attribute values
        self.ss = results['ss']
        self.shell_cmd = results['shell_cmd']
        self.varNames = results['varNames']
        self.stateNames = results['states']
        self.shockNames = results['shockNames']
        self.locations = results['locations']

    def impulse(self,T = 51,t0=1,dropFirst = 300,shocks = None,diff = True,percent = False):

        '''Initializing an instance dynarehelper.model requires values for the following varNames:

        Args:
            T:          (int) contains the parameters of the model. Default: 51
            t0:         (int) A list of strings with the names of the exogenous shocks. Default: 1
            dropFirst:  (int) A list of strings with the names of the endogenous variables. Default: 300
            shocks:     (Pandas Series,Numpy array,list) Values for the exogenous shocks. If Series,
                            index should be a subset of the exogenous shock names. If list or Numpy
                            array,must be orderred in same order as .shockNames. Default: None
            diff:       (bool) Default: Subtract the steady state from the simulated values. True
            percent:    (bool) Multiply the simulated values by 100. Default: False


        Returns:
            None

        Attributes:
            irs:    (dict) A dictionary of Pandas DataFrames containing the simulated values for 
                        each shock.

        '''

        
        # Create an array shock values for each desired impulse response calculation
        shock_cell = np.zeros((len(self.shockNames),),dtype=np.object)
        for i in range(len(self.shockNames)):
            shock_cell[i] = np.zeros([len(self.shockNames),dropFirst+T])

        # Add shock values
        try:
            for name in self.shockNames:
                shock_cell[self.locations[name]][self.locations[name],t0+dropFirst] = shocks[name]

        except:
            if shocks is not None:
                for i,item in enumerate(self.shockNames):
                    shock_cell[i][i,t0+dropFirst] = shocks[i]
            else:
                for i,item in enumerate(self.shockNames):
                    shock_cell[i][i,t0+dropFirst] = 0.01

        # Save shock data to shocks.mat for Octave
        scipy.io.savemat('shocks.mat',{'shocks':shock_cell})

        # Create simulation.m if it doesn't exist
        octave_file()

        # Prepare command for shell
        options = ['octave','simulation.m',self.filename.replace('.mod','')+'.mat',self.moduleDir]

        # Shell command
        call(options)
        
        # Import simulated data
        sim_mat = scipy.io.loadmat('simulated_data.mat')

        # Create a dictionary of Pandas DataFrames containing computed impulse responses.
        self.irs = {}
        for shock in self.shockNames:
            self.irs[shock] = pd.DataFrame()
            shock_loc = self.locations[shock]
            
            for name in self.varNames:
                
                name_loc = self.locations[name]

                if diff == True:

                    self.irs[shock][name] = sim_mat['data'][0][shock_loc][name_loc][dropFirst:] - self.ss[name]

                else:

                    self.irs[shock][name] = sim_mat['data'][0][shock_loc][name_loc][dropFirst:]

            self.irs[shock][shock] = shock_cell[shock_loc][shock_loc][dropFirst:]

        if percent==True:
            for shock in self.shockNames:
                self.irs[shock] = 100*self.irs[shock]

    def stoch_sim(self,T=51,dropFirst=300,covMat=None,seed=None,percent=False,diff=True):

        '''Initializing an instance dynarehelper.model requires values for the following varNames:

        Args:
            T:          (int) contains the parameters of the model. Default: 51
            dropFirst:  (int) A list of strings with the names of the endogenous variables. Default: 300
            covMat:     (list or numpy Array). Shock covariance matrix. Must match order of shocknames.
                            if dim is less than number of shocks, shock(s) at end of shocks will be assumed
                            to have zero variance. Default: None
            shocks:     (Pandas Series,Numpy array,list) Values for the exogenous shocks. If Series,
                            index should be a subset of the exogenous shock names. If list or Numpy
                            array,must be orderred in same order as .shockNames. Default: None
            diff:       (bool) Default: Subtract the steady state from the simulated values. True
            percent:    (bool) Multiply the simulated values by 100. Default: False


        Returns:
            None

        Attributes:
            simulated:    (Pandas DataFrame) Simulated data.

        '''

        # Set seed for Numpy random number generator
        if seed is not None and type(seed)==int:
            np.random.seed(seed)

        # Assign values to covMat if none given
        if covMat is None:
            covMat = np.eye(len(self.shockNames))

        # Simulate exogenous shocks
        eps = np.random.multivariate_normal(mean=np.zeros(len(self.shockNames)),cov=covMat,size=[dropFirst+T]).T
        
        # Save simulated shock data to .mat file for Octave
        shock_cell = np.zeros((1,),dtype=np.object)
        shock_cell[0] = eps
        scipy.io.savemat('shocks.mat',{'shocks':shock_cell})

        # Create simulation.m if it doesn't exist
        octave_file()

        # Prepare command for shell
        options = ['octave','simulation.m',self.filename.replace('.mod','')+'.mat',self.moduleDir]

        # Shell command
        call(options)
        
        # Import simulated data
        sim_mat = scipy.io.loadmat('simulated_data.mat')

        # Create Pandas DataFrames containing simulated data.
        frameDict = {}

        for name in self.varNames:
                
            name_loc = self.locations[name]

            if diff == True:

                frameDict[name] = sim_mat['data'][0][0][name_loc][dropFirst:] - self.ss[name]

            else:

                frameDict[name] = sim_mat['data'][0][0][name_loc][dropFirst:]

        for i,name in enumerate(self.shockNames):
            frameDict[name] = eps[i][dropFirst:]


        self.simulated = pd.DataFrame(frameDict,index = np.arange(T))
        if percent==True:
            self.simulated = 100*self.simulated

    def cleanup(self,allFiles=False):
        
        '''Deletes files created by Dynare++ and optionally Octave and this program during computations.
        By default deletes all files associated with filename and ending with: '.jnl', '.dump', '_f.m', 
        and '_ff.m'.

        Args:
            allFiles:   (bool) Whether to also delete simulation.m and files ending with '.mod'
                            or '.mat'. Default: False

        Returns:
            None

        Attributes:
            None

        '''

        cleanup(self.filename,allFiles)

     

def create_dynare_file(parameters,shockNames,varNames,equations,steadyStates,covMat,order,filename):

    '''Create the dynare .mod file'.

    Args:
        parameters:     (Pandas Series) contains the parameters of the model.
        shockNames:     (list) A list of strings with the names of the exogenous shocks
        varNames:       (list) A list of strings with the names of the endogenous variables.
        equations:      (list) A list with string elements that correspond to equation lines 
                            for a Dynare++ file. No semicolons.
        steadyStates:   (list,Numpy array,Pandas Series) Initial values for the Dynare++
                            nonlinear solver. Can be a Pandas Series with index equal to varNames,
                            a numpy array with values ordered to conform with varNames,or a list
                            strings with equations to be evaluated by the solver
        covMat:         (Numpy array) Shock covariance matrix
        order:          (int) Order of approximation
        filename:       (str) nome of .mod file to be created for Dynare++ computation. If string
                            doesn't end with .mod,then .mod is appended

    Returns:
        None

    '''

    # Construct the body of the file
    text = ''
    text += 'var '
    for name in varNames:
        text += name+','

    text = text[:-1]+';\n\n'

    text += 'varexo '
    for name in shockNames:
        text += name+','

    text = text[:-1]+';\n\n'

    text += 'parameters '
    for name in parameters.index:
        text += name+','

    text = text[:-1]+';\n\n'
    
    for name in parameters.index:
        text+=name+' = '+ str(parameters[name])+';\n'

    text += '\n'

    text += 'model;\n'
    for eqn in equations:
        text += eqn+';\n'
    text += 'end;\n\n'

    text += 'initval;\n'
    try:

        for name in steadyStates.index:
            text += name+' = '+str(steadyStates[name])+';\n'

    except:

        for eqn in steadyStates:
            text += eqn+';\n'

    text += 'end;\n\n'

    text += 'vcov = [\n    '
    for row in covMat:
        for element in row:
            text+=str(element)+' '
        text = text[:-1]+';\n    '
    text+='];\n\n'
    
    text+='order = '+str(int(order))+';'

    # Write file
    with open(filename,'w') as test:
        test.write(text)


def run_dynare(filename,per=None,dropFirst=None,sim=None,centralize=False,order= None,otherOptions=None):

    '''Use Dynare++ to approximate the decision rule for the nonlinear model

    Args:
        per:            (int) Number of periods to simulate. Default: 0
        dropFirst:      (int) Number of initial simulated periods to discard. Default: 0
        sim:            (int) Number of simulation periods. Default: 0
        centralize:     (bool)=False. Whether to subtract the constant term of the decision rule (not 
                            necessarily steady state for higher order approximations). Default: False
        order:          (int) Order of approximation. Overrides self.order. Default: None
        otherOptions:   (list) additional arguments to be passed to Dynare++.  Default: None

    Returns:
        Dictionary woith the following elements:

            varNames:   (list) Names of the variables
            states:     (list) Names of the state variables of the model
            shockNames: (list) Names of the shocks in the model
            ss:         (Pandas Series) computed steady state
            locations:  (list) Ordering (starting from 1) of the endogenous variables
            shell_cmd:  (str) Shell command to drun Dynare++

    '''
    
    # Construct the shell command
    options = ['dynare++']

    if per !=None:
        options.append('--per')
        options.append(str(int(per)))

    if dropFirst !=None:
        options.append('--dropFirst')
        options.append(str(int(dropFirst)))
    
    if sim != None:
        options.append('--sim')
        options.append(str(int(sim)))

    if centralize ==True:
        options.append('--centralize')
    else:
        options.append('--no-centralize')

    if order !=None:
        options.append('--order')
        options.append(str(int(order)))

    if otherOptions !=None:
        for item in otherOptions:
            options.append(item)

    options.append(filename)
    
    # Run Dynare++
    call(options)

    # Import results
    mat = scipy.io.loadmat(filename.split('.')[0]+'.mat')
    
    varNames = []
    for item in mat['dyn_vars']:
        varNames.append(item.strip())

    stateNames = []
    for item in mat['dyn_state_vars']:
        stateNames.append(item.strip())

    shockNames = []
    for item in mat['dyn_shocks']:
        shockNames.append(item.strip())

    steadyStates = mat['dyn_ss']

    locations = {}
    for name in varNames:
        locations[name] = int(mat['dyn_i_'+name]-1)

    for name in shockNames:
        locations[name] = int(mat['dyn_i_'+name]-1)

    
    ss = pd.Series()
    for i,variable_name in enumerate(varNames):
        ss[variable_name] = steadyStates[i][0]
        
    results = {
        'varNames':varNames,
        'states':stateNames,
        'shockNames':shockNames,
        'ss':ss,
        'locations':locations,
        'shell_cmd':options
    }
    
    return results

def cleanup(filename=None,allFiles=False):

    '''Deletes files created by Dynare++ and optioanlly Octave during computations. By default
    deletes all files associated with filename and ending with: '.jnl', '.dump', '_f.m', and '_ff.m'.

    Args:
        filename:   (str) filename of simulations to delete. Default: None
        allFiles:   (bool) Whether to also delete simulation.m and files ending with '.mod'
                        or '.mat'. Default: False

    Returns:
        None

    '''

    if filename!=None:
        
        try:
            remove(filename.split('.')[0]+'.dump')
        except OSError:
            pass
        try:
            remove(filename.split('.')[0]+'.jnl')
        except OSError:
            pass
        try:
            remove(filename.split('.')[0]+'_f.m')
        except OSError:
            pass
        try:
            remove(filename.split('.')[0]+'_ff.m')
        except OSError:
            pass


        if allFiles:
            try:
                remove(filename.split('.')[0]+'.mod')
            except OSError:
                pass
            try:
                remove(filename.split('.')[0]+'.mat')
            except OSError:
                pass


    else:

        dir_name = getcwd()
        files = listdir(dir_name)

        for item in files:
            if item.endswith('.jnl') or item.endswith('.dump') or item.endswith('_f.m') or item.endswith('_ff.m'):
                remove(join(dir_name,item))

        if allFiles:
            for item in files:
                if item.endswith('.mod') or item.endswith('.mat') or item.endswith('simulation.m'):
                    remove(join(dir_name,item))


def transorm_to_logs(vars_to_transform,equations,varNames,shockNames,parameters,lead_lag_max=30):

    ''' Function for converting a model in levels into a model in logs.

    Args:
        vars_to_transform:  (list) ist of variable names to be converted to log
        equations:          (list) list of equations (strings) in Dynare++ format
        variables:           (list) list of all variables in the model
        shockNames:         (list) list of all shocks in the model
        parameters:         (list) list of all parameters in the model
        lead_lag_max:       (int) maximum number of time leads/lags in dynare++ model

    Returns:
        list of modified equations
    '''

    variables = varNames.copy()
    funs = ['erf','erfc','log','exp']

    strings = shockNames + funs+ list(parameters.index)
    

    for item in strings:
        if item in vars_to_transform:
            strings.remove(item)
            variables.append(item)

    for i,eqn in enumerate(equations):
        for v in vars_to_transform:
            eqn = eqn.replace(v,'exp('+v+')')

            for v2 in variables:
                if v in v2 and v != v2:
                    eqn = eqn.replace(v2.replace(v,'exp('+v+')'),v2)
            for string in strings:
                eqn = eqn.replace(string.replace(v,'exp('+v+')'),string)
            for k in np.arange(-lead_lag_max,lead_lag_max+1,1):
                eqn = eqn.replace('exp('+v+')('+str(k)+')','exp('+v+'('+str(k)+'))')

        equations[i] = eqn
    return equations