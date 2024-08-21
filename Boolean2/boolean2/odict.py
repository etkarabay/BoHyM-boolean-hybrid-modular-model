#
# a dictionary-like class that maintains the order of insertion
# 
# based on a recipe by Igor Ghisi located at
#
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/496761
#

from collections.abc import MutableMapping as DictMixin



class odict(DictMixin):
    
    
    def __init__(self, **kwds):
        self._keys = []
        self._data = {}
        for _key, value in list(kwds.items()):
            self[key] = value
    def __iter__(self):
        for key in self._keys:
            yield key
    def __len__(self):
        return self.keys

            
    def __setitem__(self, key, value):
        if key not in self._data:
            self._keys.append(key)
        self._data[key] = value
        
    def __getitem__(self, key):
        return self._data[key]
    
    def __delitem__(self, key):
        del self._data[key]
        self._keys.remove(key)
        
    def keys(self):
        return list(self._keys)
    
    def copy(self):
        copyDict = odict()
        copyDict._data = self._data.copy()
        copyDict._keys = self._keys[:]
        return copyDict
    
    
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()