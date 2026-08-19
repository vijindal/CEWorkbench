package org.ce.scratch;

/**
 * Frozen copy of the hand-transcribed SGTE polynomials that
 * {@code LatticeStability} carried before it was migrated to parse
 * {@code inputs/unary.dat} through {@code SgteDatabase}.
 *
 * <p>Kept solely so {@link SgteCrossCheck} retains a genuinely independent
 * second implementation. Now that the façade delegates, comparing against
 * {@code LatticeStability} would compare the database to itself and pass
 * vacuously.</p>
 *
 * <p><b>Not for production use, and not to be "fixed".</b> This is a
 * historical artifact, preserved exactly as transcribed -- including its
 * {@code t < upper} range convention, which differs from SGTE's half-open
 * {@code (lower, upper]} at exact breakpoints. The cross-check reports that
 * difference separately rather than treating it as an error.</p>
 */
final class HardcodedLatticeStability {

    private HardcodedLatticeStability() {}

    public static double ghserMo(double t) {
        if (t < 2896.00) {
            return -7746.302 + 131.9197 * t - 23.56414 * t * Math.log(t)
                    - 3.443396e-3 * t * t + 0.566283e-6 * t * t * t
                    + 65812.0 / t - 0.130927e-9 * t * t * t * t;
        }
        return -30556.41 + 283.559746 * t - 42.63829 * t * Math.log(t)
                - 4849.315e30 / Math.pow(t, 9);
    }

    /** d(ghserMo)/dT. */
    public static double dGhserMoDt(double t) {
        if (t < 2896.00) {
            return 131.9197 - 23.56414 * (Math.log(t) + 1)
                    - 2.0 * 3.443396e-3 * t + 3.0 * 0.566283e-6 * t * t
                    - 65812.0 / (t * t) - 4.0 * 0.130927e-9 * t * t * t;
        }
        return 283.559746 - 42.63829 * (Math.log(t) + 1)
                + 9.0 * 4849.315e30 / Math.pow(t, 10);
    }

    public static double gliqMo(double t) {
        if (t < 2896.00) {
            return 41831.347 - 14.694912 * t + 424.519e-24 * Math.pow(t, 7) + ghserMo(t);
        }
        return 3538.963 + 271.6697 * t - 42.63829 * t * Math.log(t);
    }

    public static double dGliqMoDt(double t) {
        if (t < 2896.00) {
            return -14.694912 + 7.0 * 424.519e-24 * Math.pow(t, 6) + dGhserMoDt(t);
        }
        return 271.6697 - 42.63829 * (Math.log(t) + 1);
    }

    public static double gfccMo(double t) {
        return 15200.0 + 0.63 * t + ghserMo(t);
    }

    public static double dGfccMoDt(double t) {
        return 0.63 + dGhserMoDt(t);
    }

    public static double ghcpMo(double t) {
        return 11550.0 + ghserMo(t);
    }

    public static double dGhcpMoDt(double t) {
        return dGhserMoDt(t);
    }

    // =========================================================================
    // Nb
    // =========================================================================

    public static double ghserNb(double t) {
        if (t < 2750.00) {
            return -8519.353 + 142.045475 * t - 26.4711 * t * Math.log(t)
                    + 0.203475e-3 * t * t - 0.35012e-6 * t * t * t
                    + 93399.0 / t;
        }
        return -37669.3 + 271.720843 * t - 41.77 * t * Math.log(t)
                + 1528.238e29 / Math.pow(t, 9);
    }

    public static double dGhserNbDt(double t) {
        if (t < 2750.00) {
            return 142.045475 - 26.4711 * (Math.log(t) + 1)
                    + 2.0 * 0.203475e-3 * t - 3.0 * 0.35012e-6 * t * t
                    - 93399.0 / (t * t);
        }
        return 271.720843 - 41.77 * (Math.log(t) + 1)
                - 9.0 * 1528.238e29 / Math.pow(t, 10);
    }

    public static double gliqNb(double t) {
        if (t < 2750.00) {
            return 29781.555 - 10.816418 * t - 306.098e-25 * Math.pow(t, 7) + ghserNb(t);
        }
        return -7499.398 + 260.756148 * t - 41.77 * t * Math.log(t);
    }

    public static double dGliqNbDt(double t) {
        if (t < 2750.00) {
            return -10.816418 - 7.0 * 306.098e-25 * Math.pow(t, 6) + dGhserNbDt(t);
        }
        return 260.756148 - 41.77 * (Math.log(t) + 1);
    }

    public static double gfccNb(double t) {
        return 13500.0 + 1.7 * t + ghserNb(t);
    }

    public static double dGfccNbDt(double t) {
        return 1.7 + dGhserNbDt(t);
    }

    public static double ghcpNb(double t) {
        return 10000.0 + 2.4 * t + ghserNb(t);
    }

    public static double dGhcpNbDt(double t) {
        return 2.4 + dGhserNbDt(t);
    }

    // =========================================================================
    // Re
    // =========================================================================

    public static double ghserRe(double t) {
        if (t < 1200.00) {
            return -7695.279 + 128.421589 * t - 24.348 * t * Math.log(t)
                    - 2.53505e-3 * t * t + 0.192818e-6 * t * t * t
                    + 32915.0 / t;
        } else if (t < 2400.00) {
            return -15775.998 + 194.667426 * t - 33.586 * t * Math.log(t)
                    + 2.24565e-3 * t * t - 0.281835e-6 * t * t * t
                    + 1376270.0 / t;
        } else if (t < 3458.00) {
            return -70882.739 + 462.110749 * t - 67.956 * t * Math.log(t)
                    + 11.84945e-3 * t * t - 0.788955e-6 * t * t * t
                    + 18075200.0 / t;
        } else if (t < 5000.00) {
            return 346325.888 - 1211.371859 * t + 140.8316548 * t * Math.log(t)
                    - 33.764567e-3 * t * t + 1.053726e-6 * t * t * t
                    - 134548866.0 / t;
        }
        return -78564.296 + 346.997842 * t - 49.519 * t * Math.log(t);
    }

    public static double dGhserReDt(double t) {
        if (t < 1200.00) {
            return 128.421589 - 24.348 * (Math.log(t) + 1)
                    - 2.0 * 2.53505e-3 * t + 3.0 * 0.192818e-6 * t * t
                    - 32915.0 / (t * t);
        } else if (t < 2400.00) {
            return 194.667426 - 33.586 * (Math.log(t) + 1)
                    + 2.0 * 2.24565e-3 * t - 3.0 * 0.281835e-6 * t * t
                    - 1376270.0 / (t * t);
        } else if (t < 3458.00) {
            return 462.110749 - 67.956 * (Math.log(t) + 1)
                    + 2.0 * 11.84945e-3 * t - 3.0 * 0.788955e-6 * t * t
                    - 18075200.0 / (t * t);
        } else if (t < 5000.00) {
            return -1211.371859 + 140.8316548 * (Math.log(t) + 1)
                    - 2.0 * 33.764567e-3 * t + 3.0 * 1.053726e-6 * t * t
                    + 134548866.0 / (t * t);
        }
        return 346.997842 - 49.519 * (Math.log(t) + 1);
    }

    public static double gliqRe(double t) {
        if (t < 1200.00) {
            return 16125.604 + 122.076209 * t - 24.348 * t * Math.log(t)
                    - 2.53505e-3 * t * t + 0.192818e-6 * t * t * t
                    + 32915.0 / t;
        } else if (t < 2000.00) {
            return 8044.885 + 188.322047 * t - 33.586 * t * Math.log(t)
                    + 2.24565e-3 * t * t - 0.281835e-6 * t * t * t
                    + 1376270.0 / t;
        } else if (t < 3458.00) {
            return 568842.665 - 2527.838455 * t + 314.1788975 * t * Math.log(t)
                    - 89.39817e-3 * t * t + 3.92854e-6 * t * t * t
                    - 163100987.0 / t;
        }
        return -39044.888 + 335.723691 * t - 49.519 * t * Math.log(t);
    }

    public static double dGliqReDt(double t) {
        if (t < 1200.00) {
            return 122.076209 - 24.348 * (Math.log(t) + 1)
                    - 2.0 * 2.53505e-3 * t + 3.0 * 0.192818e-6 * t * t
                    - 32915.0 / (t * t);
        } else if (t < 2000.00) {
            return 188.322047 - 33.586 * (Math.log(t) + 1)
                    + 2.0 * 2.24565e-3 * t - 3.0 * 0.281835e-6 * t * t
                    - 1376270.0 / (t * t);
        } else if (t < 3458.00) {
            return -2527.838455 + 314.1788975 * (Math.log(t) + 1)
                    - 2.0 * 89.39817e-3 * t + 3.0 * 3.92854e-6 * t * t
                    + 163100987.0 / (t * t);
        }
        return 335.723691 - 49.519 * (Math.log(t) + 1);
    }

    public static double gbccRe(double t) {
        return 17000.0 - 3.7 * t + ghserRe(t);
    }

    public static double dGbccReDt(double t) {
        return -3.7 + dGhserReDt(t);
    }

    public static double gfccRe(double t) {
        return 11000.0 - 1.5 * t + ghserRe(t);
    }

    public static double dGfccReDt(double t) {
        return -1.5 + dGhserReDt(t);
    }

    // =========================================================================
    // Ta
    // =========================================================================

    public static double ghserTa(double t) {
        if (t < 1300.00) {
            return -7285.889 + 119.139857 * t - 23.7592624 * t * Math.log(t)
                    - 2.623033e-3 * t * t + 0.170109e-6 * t * t * t
                    - 3293.0 / t;
        } else if (t < 2500.00) {
            return -22389.955 + 243.88676 * t - 41.137088 * t * Math.log(t)
                    + 6.167572e-3 * t * t - 0.655136e-6 * t * t * t
                    + 2429586.0 / t;
        } else if (t < 3290.00) {
            return 229382.886 - 722.59722 * t + 78.5244752 * t * Math.log(t)
                    - 17.983376e-3 * t * t + 0.195033e-6 * t * t * t
                    - 93813648.0 / t;
        }
        return -1042384.014 + 2985.491246 * t - 362.1591318 * t * Math.log(t)
                + 43.117795e-3 * t * t - 1.055148e-6 * t * t * t
                + 554714342.0 / t;
    }

    public static double dGhserTaDt(double t) {
        if (t < 1300.00) {
            return 119.139857 - 23.7592624 * (Math.log(t) + 1)
                    - 2.0 * 2.623033e-3 * t + 3.0 * 0.170109e-6 * t * t
                    + 3293.0 / (t * t);
        } else if (t < 2500.00) {
            return 243.88676 - 41.137088 * (Math.log(t) + 1)
                    + 2.0 * 6.167572e-3 * t - 3.0 * 0.655136e-6 * t * t
                    - 2429586.0 / (t * t);
        } else if (t < 3290.00) {
            return -722.59722 + 78.5244752 * (Math.log(t) + 1)
                    - 2.0 * 17.983376e-3 * t + 3.0 * 0.195033e-6 * t * t
                    + 93813648.0 / (t * t);
        }
        return 2985.491246 - 362.1591318 * (Math.log(t) + 1)
                + 2.0 * 43.117795e-3 * t - 3.0 * 1.055148e-6 * t * t
                - 554714342.0 / (t * t);
    }

    public static double gliqTa(double t) {
        if (t < 1000.00) {
            return 21875.086 + 111.561128 * t - 23.7592624 * t * Math.log(t)
                    - 2.623033e-3 * t * t + 0.170109e-6 * t * t * t
                    - 3293.0 / t;
        } else if (t < 3290.00) {
            return 43884.339 - 61.981795 * t + 0.0279523 * t * Math.log(t)
                    - 12.330066e-3 * t * t + 0.614599e-6 * t * t * t
                    - 3523338.0 / t;
        }
        return -6314.543 + 258.110873 * t - 41.84 * t * Math.log(t);
    }

    public static double dGliqTaDt(double t) {
        if (t < 1000.00) {
            return 111.561128 - 23.7592624 * (Math.log(t) + 1)
                    - 2.0 * 2.623033e-3 * t + 3.0 * 0.170109e-6 * t * t
                    + 3293.0 / (t * t);
        } else if (t < 3290.00) {
            return -61.981795 + 0.0279523 * (Math.log(t) + 1)
                    - 2.0 * 12.330066e-3 * t + 3.0 * 0.614599e-6 * t * t
                    + 3523338.0 / (t * t);
        }
        return 258.110873 - 41.84 * (Math.log(t) + 1);
    }

    public static double gfccTa(double t) {
        return 16000.0 + 1.7 * t + ghserTa(t);
    }

    public static double dGfccTaDt(double t) {
        return 1.7 + dGhserTaDt(t);
    }

    public static double ghcpTa(double t) {
        return 12000.0 + 2.4 * t + ghserTa(t);
    }

    public static double dGhcpTaDt(double t) {
        return 2.4 + dGhserTaDt(t);
    }

    // =========================================================================
    // Ti
    // =========================================================================

    public static double ghserTi(double t) {
        if (t < 900.00) {
            return -8059.921 + 133.615208 * t - 23.9933 * t * Math.log(t)
                    - 4.777975e-3 * t * t + 0.106716e-6 * t * t * t
                    + 72636.0 / t;
        } else if (t < 1155.00) {
            return -7811.815 + 132.988068 * t - 23.9887 * t * Math.log(t)
                    - 4.2033e-3 * t * t - 0.090876e-6 * t * t * t
                    + 42680.0 / t;
        } else if (t < 1941.00) {
            return 908.837 + 66.976538 * t - 14.9466 * t * Math.log(t)
                    - 8.1465e-3 * t * t + 0.202715e-6 * t * t * t
                    - 1477660.0 / t;
        }
        return -124526.786 + 638.806871 * t - 87.2182461 * t * Math.log(t)
                + 8.204849e-3 * t * t - 0.304747e-6 * t * t * t
                + 36699805.0 / t;
    }

    public static double dGhserTiDt(double t) {
        if (t < 900.00) {
            return 133.615208 - 23.9933 * (Math.log(t) + 1)
                    - 2.0 * 4.777975e-3 * t + 3.0 * 0.106716e-6 * t * t
                    - 72636.0 / (t * t);
        } else if (t < 1155.00) {
            return 132.988068 - 23.9887 * (Math.log(t) + 1)
                    - 2.0 * 4.2033e-3 * t - 3.0 * 0.090876e-6 * t * t
                    - 42680.0 / (t * t);
        } else if (t < 1941.00) {
            return 66.976538 - 14.9466 * (Math.log(t) + 1)
                    - 2.0 * 8.1465e-3 * t + 3.0 * 0.202715e-6 * t * t
                    + 1477660.0 / (t * t);
        }
        return 638.806871 - 87.2182461 * (Math.log(t) + 1)
                + 2.0 * 8.204849e-3 * t - 3.0 * 0.304747e-6 * t * t
                - 36699805.0 / (t * t);
    }

    public static double gbccTi(double t) {
        if (t < 1155.00) {
            return -1272.064 + 134.71418 * t - 25.5768 * t * Math.log(t)
                    - 0.663845e-3 * t * t - 0.278803e-6 * t * t * t
                    + 7208.0 / t;
        } else if (t < 1941.00) {
            return 6667.385 + 105.366379 * t - 22.3771 * t * Math.log(t)
                    + 1.21707e-3 * t * t - 0.84534e-6 * t * t * t
                    - 2002750.0 / t;
        }
        // Reference has no branch beyond 4000 K for GBCCTI (Module leaves it
        // unassigned); mirror that boundary rather than silently extrapolate.
        if (t < 4000.00) {
            return 26483.26 - 182.426471 * t + 19.0900905 * t * Math.log(t)
                    - 22.00832e-3 * t * t + 1.228863e-6 * t * t * t
                    + 1400501.0 / t;
        }
        throw new IllegalArgumentException("gbccTi undefined above 4000 K (T=" + t + ")");
    }

    public static double dGbccTiDt(double t) {
        if (t < 1155.00) {
            return 134.71418 - 25.5768 * (Math.log(t) + 1)
                    - 2.0 * 0.663845e-3 * t - 3.0 * 0.278803e-6 * t * t
                    - 7208.0 / (t * t);
        } else if (t < 1941.00) {
            return 105.366379 - 22.3771 * (Math.log(t) + 1)
                    + 2.0 * 1.21707e-3 * t - 3.0 * 0.84534e-6 * t * t
                    + 2002750.0 / (t * t);
        }
        if (t < 4000.00) {
            return -182.426471 + 19.0900905 * (Math.log(t) + 1)
                    - 2.0 * 22.00832e-3 * t + 3.0 * 1.228863e-6 * t * t
                    - 1400501.0 / (t * t);
        }
        throw new IllegalArgumentException("dGbccTiDt undefined above 4000 K (T=" + t + ")");
    }

    public static double gliqTi(double t) {
        if (t < 900.00) {
            return 4134.494 + 126.63427 * t - 23.9933 * t * Math.log(t)
                    - 4.777975e-3 * t * t + 0.106716e-6 * t * t * t
                    + 72636.0 / t;
        } else if (t < 1155.00) {
            return 4382.601 + 126.00713 * t - 23.9887 * t * Math.log(t)
                    - 4.2033e-3 * t * t - 0.090876e-6 * t * t * t
                    + 42680.0 / t;
        } else if (t < 1300.00) {
            return 13103.253 + 59.9956 * t - 14.9466 * t * Math.log(t)
                    - 8.1465e-3 * t * t + 0.202715e-6 * t * t * t
                    - 1477660.0 / t;
        } else if (t < 1941.00) {
            return 369519.198 - 2554.0225 * t + 342.059267 * t * Math.log(t)
                    - 163.409355e-3 * t * t + 12.457117e-6 * t * t * t
                    - 67034516.0 / t;
        }
        return -19887.066 + 298.7367 * t - 46.29 * t * Math.log(t);
    }

    public static double dGliqTiDt(double t) {
        if (t < 900.00) {
            return 126.63427 - 23.9933 * (Math.log(t) + 1)
                    - 2.0 * 4.777975e-3 * t + 3.0 * 0.106716e-6 * t * t
                    - 72636.0 / (t * t);
        } else if (t < 1155.00) {
            return 126.00713 - 23.9887 * (Math.log(t) + 1)
                    - 2.0 * 4.2033e-3 * t - 3.0 * 0.090876e-6 * t * t
                    - 42680.0 / (t * t);
        } else if (t < 1300.00) {
            return 59.9956 - 14.9466 * (Math.log(t) + 1)
                    - 2.0 * 8.1465e-3 * t + 3.0 * 0.202715e-6 * t * t
                    + 1477660.0 / (t * t);
        } else if (t < 1941.00) {
            return -2554.0225 + 342.059267 * (Math.log(t) + 1)
                    - 2.0 * 163.409355e-3 * t + 3.0 * 12.457117e-6 * t * t
                    + 67034516.0 / (t * t);
        }
        return 298.7367 - 46.29 * (Math.log(t) + 1);
    }

    public static double gfccTi(double t) {
        return 6000.0 - 0.1 * t + ghserTi(t);
    }

    public static double dGfccTiDt(double t) {
        return -0.1 + dGhserTiDt(t);
    }

    // =========================================================================
    // V
    // =========================================================================

    public static double ghserV(double t) {
        if (t < 790.00) {
            return -7930.43 + 133.346053 * t - 24.134 * t * Math.log(t)
                    - 3.098e-3 * t * t + 0.12175e-6 * t * t * t
                    + 69460.0 / t;
        } else if (t < 2183.00) {
            return -7967.842 + 143.291093 * t - 25.9 * t * Math.log(t)
                    + 0.0625e-3 * t * t - 0.68e-6 * t * t * t;
        }
        // Reference has no branch beyond 4000 K for GHSERV (Module leaves it
        // unassigned); mirror that boundary rather than silently extrapolate.
        if (t < 4000.00) {
            return -41689.864 + 321.140783 * t - 47.43 * t * Math.log(t)
                    + 644.389e29 / Math.pow(t, 9);
        }
        throw new IllegalArgumentException("ghserV undefined above 4000 K (T=" + t + ")");
    }

    public static double dGhserVDt(double t) {
        if (t < 790.00) {
            return 133.346053 - 24.134 * (Math.log(t) + 1)
                    - 2.0 * 3.098e-3 * t + 3.0 * 0.12175e-6 * t * t
                    - 69460.0 / (t * t);
        } else if (t < 2183.00) {
            return 143.291093 - 25.9 * (Math.log(t) + 1)
                    + 2.0 * 0.0625e-3 * t - 3.0 * 0.68e-6 * t * t;
        }
        if (t < 4000.00) {
            return 321.140783 - 47.43 * (Math.log(t) + 1)
                    - 9.0 * 644.389e29 / Math.pow(t, 10);
        }
        throw new IllegalArgumentException("dGhserVDt undefined above 4000 K (T=" + t + ")");
    }

    public static double ghcpV(double t) {
        if (t < 790.00) {
            return -3930.43 + 135.746053 * t - 24.134 * t * Math.log(t)
                    - 3.098e-3 * t * t + 0.12175e-6 * t * t * t
                    + 69460.0 / t;
        } else if (t < 2183.00) {
            return -3967.842 + 145.691093 * t - 25.9 * t * Math.log(t)
                    + 0.0625e-3 * t * t - 0.68e-6 * t * t * t;
        }
        if (t < 4000.00) {
            return -37689.864 + 323.540783 * t - 47.43 * t * Math.log(t)
                    + 644.389e29 / Math.pow(t, 9);
        }
        throw new IllegalArgumentException("ghcpV undefined above 4000 K (T=" + t + ")");
    }

    public static double dGhcpVDt(double t) {
        if (t < 790.00) {
            return 135.746053 - 24.134 * (Math.log(t) + 1)
                    - 2.0 * 3.098e-3 * t + 3.0 * 0.12175e-6 * t * t
                    - 69460.0 / (t * t);
        } else if (t < 2183.00) {
            return 145.691093 - 25.9 * (Math.log(t) + 1)
                    + 2.0 * 0.0625e-3 * t - 3.0 * 0.68e-6 * t * t;
        }
        if (t < 4000.00) {
            return 323.540783 - 47.43 * (Math.log(t) + 1)
                    - 9.0 * 644.389e29 / Math.pow(t, 10);
        }
        throw new IllegalArgumentException("dGhcpVDt undefined above 4000 K (T=" + t + ")");
    }

    public static double gliqV(double t) {
        if (t < 790.00) {
            return 12833.687 + 123.890501 * t - 24.134 * t * Math.log(t)
                    - 3.098e-3 * t * t + 0.12175e-6 * t * t * t
                    + 69460.0 / t - 519.136e-24 * Math.pow(t, 7);
        } else if (t < 2183.00) {
            return 12796.275 + 133.835541 * t - 25.9 * t * Math.log(t)
                    + 0.0625e-3 * t * t - 0.68e-6 * t * t * t
                    - 519.136e-24 * Math.pow(t, 7);
        }
        if (t < 4000.00) {
            return -19617.51 + 311.055983 * t - 47.43 * t * Math.log(t);
        }
        throw new IllegalArgumentException("gliqV undefined above 4000 K (T=" + t + ")");
    }

    public static double dGliqVDt(double t) {
        if (t < 790.00) {
            return 123.890501 - 24.134 * (Math.log(t) + 1)
                    - 2.0 * 3.098e-3 * t + 3.0 * 0.12175e-6 * t * t
                    - 69460.0 / (t * t) - 7.0 * 519.136e-24 * Math.pow(t, 6);
        } else if (t < 2183.00) {
            return 133.835541 - 25.9 * (Math.log(t) + 1)
                    + 2.0 * 0.0625e-3 * t - 3.0 * 0.68e-6 * t * t
                    - 7.0 * 519.136e-24 * Math.pow(t, 6);
        }
        if (t < 4000.00) {
            return 311.055983 - 47.43 * (Math.log(t) + 1);
        }
        throw new IllegalArgumentException("dGliqVDt undefined above 4000 K (T=" + t + ")");
    }

    // =========================================================================
    // W
    // =========================================================================

    public static double ghserW(double t) {
        if (t < 3695.00) {
            return -7646.311 + 130.4 * t - 24.1 * t * Math.log(t)
                    - 1.936e-3 * t * t + 0.207e-6 * t * t * t
                    + 44500.0 / t - 0.0533e-9 * t * t * t * t;
        }
        return -82868.801 + 389.362335 * t - 54.0 * t * Math.log(t)
                + 1528.621e30 / Math.pow(t, 9);
    }

    public static double dGhserWDt(double t) {
        if (t < 3695.00) {
            return 130.4 - 24.1 * (Math.log(t) + 1)
                    - 2.0 * 1.936e-3 * t + 3.0 * 0.207e-6 * t * t
                    - 44500.0 / (t * t) - 4.0 * 0.0533e-9 * t * t * t;
        }
        return 389.362335 - 54.0 * (Math.log(t) + 1)
                - 9.0 * 1528.621e30 / Math.pow(t, 10);
    }

    public static double gliqW(double t) {
        if (t < 3695.00) {
            return 52160.584 - 14.10999 * t - 2713.468e-27 * Math.pow(t, 7) + ghserW(t);
        }
        return -30436.051 + 375.175 * t - 54.0 * t * Math.log(t);
    }

    public static double dGliqWDt(double t) {
        if (t < 3695.00) {
            return -14.10999 - 7.0 * 2713.468e-27 * Math.pow(t, 6) + dGhserWDt(t);
        }
        return 375.175 - 54.0 * (Math.log(t) + 1);
    }

    public static double gfccW(double t) {
        return 19300.0 + 0.63 * t + ghserW(t);
    }

    public static double dGfccWDt(double t) {
        return 0.63 + dGhserWDt(t);
    }

    public static double ghcpW(double t) {
        return 14750.0 + ghserW(t);
    }

    public static double dGhcpWDt(double t) {
        return dGhserWDt(t);
    }

    // =========================================================================
    // Zr
    // =========================================================================

    public static double ghserZr(double t) {
        if (t < 2128.00) {
            return -7827.595 + 125.64905 * t - 24.1618 * t * Math.log(t)
                    - 4.37791e-3 * t * t + 34971.0 / t;
        }
        return -26085.921 + 262.724183 * t - 42.144 * t * Math.log(t)
                - 1342.896e28 / Math.pow(t, 9);
    }

    public static double dGhserZrDt(double t) {
        if (t < 2128.00) {
            return 125.64905 - 24.1618 * (Math.log(t) + 1)
                    - 2.0 * 4.37791e-3 * t - 34971.0 / (t * t);
        }
        return 262.724183 - 42.144 * (Math.log(t) + 1)
                + 9.0 * 1342.896e28 / Math.pow(t, 10);
    }

    public static double gbccZr(double t) {
        if (t < 2128.00) {
            return -525.539 + 124.9457 * t - 25.607406 * t * Math.log(t)
                    - 0.340084e-3 * t * t - 0.009729e-6 * t * t * t
                    + 25233.0 / t - 0.076143e-9 * t * t * t * t;
        }
        return -30705.955 + 264.284163 * t - 42.144 * t * Math.log(t)
                + 1276.058e29 / Math.pow(t, 9);
    }

    public static double dGbccZrDt(double t) {
        if (t < 2128.00) {
            return 124.9457 - 25.607406 * (Math.log(t) + 1)
                    - 2.0 * 0.340084e-3 * t - 3.0 * 0.009729e-6 * t * t
                    - 25233.0 / (t * t) - 4.0 * 0.076143e-9 * t * t * t;
        }
        return 264.284163 - 42.144 * (Math.log(t) + 1)
                - 9.0 * 1276.058e29 / Math.pow(t, 10);
    }

    public static double gliqZr(double t) {
        if (t < 2128.00) {
            return 10320.095 + 116.568238 * t - 24.1618 * t * Math.log(t)
                    - 4.37791e-3 * t * t + 34971.0 / t
                    + 1627.5e-25 * Math.pow(t, 7);
        }
        return -8281.26 + 253.812609 * t - 42.144 * t * Math.log(t);
    }

    public static double dGliqZrDt(double t) {
        if (t < 2128.00) {
            return 116.568238 - 24.1618 * (Math.log(t) + 1)
                    - 2.0 * 4.37791e-3 * t - 34971.0 / (t * t)
                    + 7.0 * 1627.5e-25 * Math.pow(t, 6);
        }
        return 253.812609 - 42.144 * (Math.log(t) + 1);
    }

    // =========================================================================
    // Dispatch: calG0 equivalent
    // =========================================================================

    /**
     * Pure-element reference energy for one element in one structure at
     * temperature {@code t} (K) -- port of the reference's {@code calG0}.
     *
     * @throws IllegalArgumentException if the (element, phase) pair is not
     *         defined in the reference data (e.g. vanadium has no FCC_A1
     *         entry) -- never silently returns zero for a gap in coverage.
     */
    public static double g0(String element, String phase, double t) {
        return switch (element) {
            case "Mo" -> switch (phase) {
                case "LIQUID" -> gliqMo(t);
                case "BCC_A2" -> ghserMo(t);
                case "FCC_A1" -> gfccMo(t);
                case "HCP_A3" -> ghcpMo(t);
                default -> throw undefined(element, phase);
            };
            case "Nb" -> switch (phase) {
                case "LIQUID" -> gliqNb(t);
                case "BCC_A2" -> ghserNb(t);
                case "FCC_A1" -> gfccNb(t);
                case "HCP_A3" -> ghcpNb(t);
                default -> throw undefined(element, phase);
            };
            case "Re" -> switch (phase) {
                case "LIQUID" -> gliqRe(t);
                case "BCC_A2" -> gbccRe(t);
                case "FCC_A1" -> gfccRe(t);
                case "HCP_A3" -> ghserRe(t);
                default -> throw undefined(element, phase);
            };
            case "Ta" -> switch (phase) {
                case "LIQUID" -> gliqTa(t);
                case "BCC_A2" -> ghserTa(t);
                case "FCC_A1" -> gfccTa(t);
                case "HCP_A3" -> ghcpTa(t);
                default -> throw undefined(element, phase);
            };
            case "Ti" -> switch (phase) {
                case "LIQUID" -> gliqTi(t);
                case "BCC_A2" -> gbccTi(t);
                case "FCC_A1" -> gfccTi(t);
                case "HCP_A3" -> ghserTi(t);
                default -> throw undefined(element, phase);
            };
            case "V" -> switch (phase) {
                case "LIQUID" -> gliqV(t);
                case "BCC_A2" -> ghserV(t);
                case "HCP_A3" -> ghcpV(t);
                default -> throw undefined(element, phase);
            };
            case "W" -> switch (phase) {
                case "LIQUID" -> gliqW(t);
                case "BCC_A2" -> ghserW(t);
                case "FCC_A1" -> gfccW(t);
                case "HCP_A3" -> ghcpW(t);
                default -> throw undefined(element, phase);
            };
            case "Zr" -> switch (phase) {
                case "LIQUID" -> gliqZr(t);
                case "BCC_A2" -> gbccZr(t);
                case "HCP_A3" -> ghserZr(t);
                default -> throw undefined(element, phase);
            };
            default -> throw undefined(element, phase);
        };
    }

    /**
     * Analytic {@code d(g0)/dT} for one element in one structure -- mirrors
     * {@link #g0}'s dispatch exactly, routing to the matching {@code dXxxDt}
     * derivative of whichever {@code Xxx} function {@link #g0} would have
     * called at the same {@code (element, phase)}.
     *
     * @throws IllegalArgumentException under the same conditions as {@link #g0}.
     */
    public static double dG0Dt(String element, String phase, double t) {
        return switch (element) {
            case "Mo" -> switch (phase) {
                case "LIQUID" -> dGliqMoDt(t);
                case "BCC_A2" -> dGhserMoDt(t);
                case "FCC_A1" -> dGfccMoDt(t);
                case "HCP_A3" -> dGhcpMoDt(t);
                default -> throw undefined(element, phase);
            };
            case "Nb" -> switch (phase) {
                case "LIQUID" -> dGliqNbDt(t);
                case "BCC_A2" -> dGhserNbDt(t);
                case "FCC_A1" -> dGfccNbDt(t);
                case "HCP_A3" -> dGhcpNbDt(t);
                default -> throw undefined(element, phase);
            };
            case "Re" -> switch (phase) {
                case "LIQUID" -> dGliqReDt(t);
                case "BCC_A2" -> dGbccReDt(t);
                case "FCC_A1" -> dGfccReDt(t);
                case "HCP_A3" -> dGhserReDt(t);
                default -> throw undefined(element, phase);
            };
            case "Ta" -> switch (phase) {
                case "LIQUID" -> dGliqTaDt(t);
                case "BCC_A2" -> dGhserTaDt(t);
                case "FCC_A1" -> dGfccTaDt(t);
                case "HCP_A3" -> dGhcpTaDt(t);
                default -> throw undefined(element, phase);
            };
            case "Ti" -> switch (phase) {
                case "LIQUID" -> dGliqTiDt(t);
                case "BCC_A2" -> dGbccTiDt(t);
                case "FCC_A1" -> dGfccTiDt(t);
                case "HCP_A3" -> dGhserTiDt(t);
                default -> throw undefined(element, phase);
            };
            case "V" -> switch (phase) {
                case "LIQUID" -> dGliqVDt(t);
                case "BCC_A2" -> dGhserVDt(t);
                case "HCP_A3" -> dGhcpVDt(t);
                default -> throw undefined(element, phase);
            };
            case "W" -> switch (phase) {
                case "LIQUID" -> dGliqWDt(t);
                case "BCC_A2" -> dGhserWDt(t);
                case "FCC_A1" -> dGfccWDt(t);
                case "HCP_A3" -> dGhcpWDt(t);
                default -> throw undefined(element, phase);
            };
            case "Zr" -> switch (phase) {
                case "LIQUID" -> dGliqZrDt(t);
                case "BCC_A2" -> dGbccZrDt(t);
                case "HCP_A3" -> dGhserZrDt(t);
                default -> throw undefined(element, phase);
            };
            default -> throw undefined(element, phase);
        };
    }

    private static IllegalArgumentException undefined(String element, String phase) {
        return new IllegalArgumentException(
                "No lattice-stability data for element=" + element + ", phase=" + phase);
    }
}
