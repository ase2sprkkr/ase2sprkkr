def fortran_format(n, format='{:.12e}'):
    la = ('{'+format+'}').format(float(n))
    a = la.lstrip()
    leading = la[:len(la)-len(a)]
    la = a
    a = la.rstrip()
    trailing = la[len(a):]
    if 'E' in format:
        e = a.find('E')
    elif 'e' in format:
        e = a.find('e')
    else:
        raise("No E in fortran format string: " + format)
    return leading + '0.{}{}{}{:02d}'.format(a[0],a[2:e],a[e:e+2],abs(int(a[e+1:])*1+1)) + trailing
