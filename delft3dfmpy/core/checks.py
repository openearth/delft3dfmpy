import logging

logger = logging.getLogger(__name__)


def check_argument(argument, name, types, columns=None):
    if not isinstance(argument, types):
        raise TypeError(f'Expected argument type {types} for variable {name}, got {type(argument)}.')

    # Check columns
    if columns is not None:
        if isinstance(columns, str):
            columns = [columns]
        missing = ', '.join([col for col in columns if col not in argument.columns])
        if any(missing):
            raise KeyError(f'Missing column(s) "{missing}"')

def warn_dataframe_not_empty(dataframe):
    if not dataframe.empty:
        logger.warning('The dataframe is not empty. The data are overwritten.')

def check_dictionary(dct, required, choice=[]):
    if required is not None:
        if isinstance(required, str):
            required = [required]
        for key in required:
            if not key in dct.keys():
                raise KeyError(f'Key "{key}" missing from dictionary.')

    if any(choice):
        for i, option in enumerate(choice):
            if not isinstance(option, (list, set)):
                choice[i] = [option]
        
        if not any([all([key in dct.keys() for key in option]) for option in choice]):
            options = [', '.join(option) for option in choice]
            raise KeyError('One of the following key sets should be present in the dictionary:\n - '.format('\n - '.join(options)))
            
