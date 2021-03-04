
def cart2Radec(cart):

    """
    Extension to Cart object to correctly convert theta and phi to Right Ascension and Declination
    inputs:
    cart [Vector3D Object] - cart object to convert
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

def aries2greenwich(cart, JD):

    """
    Rotates the coordinate system of cart from the First Point of Aries to the location of 
    Greenwich at a given julian date
    inputs:
    cart [Vector3D Obj] - vector to rotate
    JD [float] - Julian date (to find the position of Greenwich)
    """

    # angle is negative in the notes
    shift = Angle(360 - jd2Angle(JD).deg, deg=True)

    return cart.rotate(shift, axis="z")
