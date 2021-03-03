
def cart2Radec(cart):

    """
    Extension to Cart object to correctly convert theta and phi to Right Ascension and Declination
    inputs:
    cart [Cart Object] - cart object to convert
    returns:
    ra [Right Ascension Obj] - Right ascension of cart
    dec [float] - Declination of cart
    """
    if not hasattr(cart, "theta") or not hasattr(cart, "phi"):
        print("[WARNING] Cartesian object does not have angles!")
        return None, None

    dec = 90 - cart.theta.deg
    ra  = cart.phi.deg

    return ra, dec