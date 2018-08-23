#
# Stolen from here: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
#
def progressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()


class InvalidCoordinatesError(Exception):

    def __init__(self, message, errors = None):
        super().__init__(message)

        self.errors = errors

def checkCoordinates(lat,lon):

    if lat > 90 or lat < -90:
        raise InvalidCoordinatesError("latitude " + str(lat) + " not allowed")
    
    if lon > 180 or lon < -180:
        raise InvalidCoordinatesError("longitude " + str(lon) + " not allowed")
    
    return True