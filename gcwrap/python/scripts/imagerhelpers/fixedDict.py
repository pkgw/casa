
class fixedDict(dict):
    """
        This allows for imposing a set keys. When using as parameters it forces 
        the user to use only allowed keys
    """
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
       

    def __delitem__(self, key):
        dict.__delitem__(self, key)

    def __setitem__(self, key, item):
        if key not in list(self.keys()): 
            raise KeyError("The key {} is not defined.".format(key))
        dict.__setitem__(self, key, item)

    def __repr__(self):
        dictrepr = dict.__repr__(self)
        return dictrepr
#        return '%s(%s)' % (type(self).__name__, dictrepr)
         
    def clear(self):
        dict.clear(self)

    def copy(self):
        copydict = dict.copy(self)
        return fixedDict(copydict)
    
    def items(self):
        return list(zip(list(self.keys()), list(self.values())))

    def keys(self):
        return dict.keys(self)

    def update(self, *args, **kwargs):
        for key, value  in dict(*args, **kwargs).items():
            dict.__setitem__(self, key,value)

   
