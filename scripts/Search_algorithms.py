

def JumpSearch (lys, val):

    
    '''This function returns the position of the element in the list lys that contains the string pattern lys. If no match 
    
    Usage: lys = list of strings to search through; val = string pattern to search for
    
    Warning: This function only works when the beginning of the string matches val'''
    
    import math
    import re
    
    length = len(lys)
    jump = int(math.sqrt(length))
    left, right = 0, 0
    index_list = sorted([lys[left],val,lys[right]])
    p = re.compile(re.escape(val))
    while left < length and index_list.index(lys[left]) <= index_list.index(val):
        right = min(length - 1, left + jump)
        index_list = sorted([lys[left],val,lys[right]])
        if index_list.index(lys[left]) <= index_list.index(val) and index_list.index(lys[right]) >= index_list.index(val):
            break
        left += jump;
        
    if left >= length or index_list.index(lys[left]) > index_list.index(val):
        return -1
    right = min(length-1, right)
    i = left
    #index_list = sorted([lys[i],val])
    while i <= right:
        index_list = sorted([lys[i],val])
        #print(p.search(lys[i]), lys[i])
        if p.match(lys[i]):
            return i
        i += 1
      
    return -1


def BinarySearch(lys, val):
    
    '''requires re'''
    
    '''This function returns the position of the element in the list lys that contains the string pattern lys. If no match 
    
    Usage: lys = list of strings to search through; val = string pattern to search for
    
    Warning: This function only works when the beginning of the string matches val'''
    
    import re
        
    first = 0
    last = len(lys)-1
    index = -1
    
    
    p = re.compile(re.escape(val))
    
    
    while (first <= last) and (index == -1):
        mid = round((first+last)/2)
        
        index_list = sorted([lys[mid],val])
        
        if p.match(lys[mid]):
            index = mid
        else:
            
            if index_list.index(val)<index_list.index(lys[mid]):
                last = mid -1
                
            else:
                first = mid +1
                
    return index