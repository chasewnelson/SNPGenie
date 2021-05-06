###############################################################################
# FUNCTION to calculate Tajima's D and its significance level, based on
# Tajima 1989 Statistical Method for Testing... Genetics 123:585-595

# AUTHOR: Chase W. Nelson
# Copyright (C) 2016 Chase W. Nelson
# DATE CREATED: April 2016

# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com

# AFFILIATION: Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA

# CITATION1: SNPGenie, https://github.com/chasewnelson/snpgenie
# CITATION2: Nelson CW, Moncla LH, Hughes AL (2015) SNPGenie: estimating evolutionary 
#	  parameters to detect natural selection using pooled next-generation sequencing data. 
#	  Bioinformatics 31(22):3709-11, doi: 10.1093/bioinformatics/btv449.

# INPUT TO FUNCTION:
# n in number of sequences (median coverage?)
# S is the number of segregating sites (not corrected)
# pi is the nucleotide diversity of the sequence

###############################################################################
# EXAMPLE:
#Tajima.D(1000, 15, 0.05)

# RETURNS:
#D  sig.level
#-2.142552	0.01

###############################################################################
Tajima.D <- function(n, S, pi) { 
  a1 <- 0
  a2 <- 0
  
  for(i in 1:(n-1)) {
    a1 <- a1 + (1/i)
    a2 <- a2 + (1/(i^2))
  }
  
  b1 <- (n + 1) / (3*(n - 1))
  b2 <- 2*(n^2 + n + 3) / (9*n*(n - 1))
  
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n + 2)/(a1 * n)) + (a2 / a1^2)
  
  e1 <- c1 / a1;
  e2 <- c2 / (a1^2 + a2)
  
  V <- (e1 * S) + (e2 * S * (S - 1))
  std.err <- sqrt(V)
  
  D <- (pi - (S / a1)) / std.err

  sig.level <- 'n.s.'
  
  if(n >= 1000) {
    if(D <= -2.369 | D >= 3.722) {
      sig.level <- 0.001
    } else if(D <= -2.062 | D >= 2.887) {
      sig.level <- 0.01
    } else if(D <= -1.715 | D >= 2.150) {
      sig.level <- 0.05
    } else if(D <= -1.505 | D >= 1.772) {
      sig.level <- 0.1
    }
  } else if(n >= 800) {
    if(D <= -2.382 | D >= 3.694) {
      sig.level <- 0.001
    } else if(D <= -2.072 | D >= 2.873) {
      sig.level <- 0.01
    } else if(D <= -1.721 | D >= 2.143) {
      sig.level <- 0.05
    } else if(D <= -1.510 | D >= 1.769) {
      sig.level <- 0.1
    }
  } else if(n >= 600) {
    if(D <= -2.398 | D >= 3.657) {
      sig.level <- 0.001
    } else if(D <= -2.084 | D >= 2.853) {
      sig.level <- 0.01
    } else if(D <= -1.728 | D >= 2.135) {
      sig.level <- 0.05
    } else if(D <= -1.515 | D >= 1.765) {
      sig.level <- 0.1
    }
  } else if(n >= 500) {
    if(D <= -2.409 | D >= 3.632) {
      sig.level <- 0.001
    } else if(D <= -2.092 | D >= 2.840) {
      sig.level <- 0.01
    } else if(D <= -1.734 | D >= 2.130) {
      sig.level <- 0.05
    } else if(D <= -1.519 | D >= 1.763) {
      sig.level <- 0.1
    }
  } else if(n >= 450) {
    if(D <= -2.415 | D >= 3.617) {
      sig.level <- 0.001
    } else if(D <= -2.096 | D >= 2.833) {
      sig.level <- 0.01
    } else if(D <= -1.737 | D >= 2.127) {
      sig.level <- 0.05
    } else if(D <= -1.521 | D >= 1.761) {
      sig.level <- 0.1
    }
  } else if(n >= 400) {
    if(D <= -2.422 | D >= 3.600) {
      sig.level <- 0.001
    } else if(D <= -2.101 | D >= 2.824) {
      sig.level <- 0.01
    } else if(D <= -1.740 | D >= 2.123) {
      sig.level <- 0.05
    } else if(D <= -1.523 | D >= 1.759) {
      sig.level <- 0.1
    }
  } else if(n >= 350) {
    if(D <= -2.430 | D >= 3.581) {
      sig.level <- 0.001
    } else if(D <= -2.107 | D >= 2.814) {
      sig.level <- 0.01
    } else if(D <= -1.744 | D >= 2.119) {
      sig.level <- 0.05
    } else if(D <= -1.526 | D >= 1.757) {
      sig.level <- 0.1
    }
  } else if(n >= 300) {
    if(D <= -2.439 | D >= 3.558) {
      sig.level <- 0.001
    } else if(D <= -2.114 | D >= 2.802) {
      sig.level <- 0.01
    } else if(D <= -1.748 | D >= 2.114) {
      sig.level <- 0.05
    } else if(D <= -1.530 | D >= 1.755) {
      sig.level <- 0.1
    }
  } else if(n >= 250) {
    if(D <= -2.449 | D >= 3.529) {
      sig.level <- 0.001
    } else if(D <= -2.122 | D >= 2.787) {
      sig.level <- 0.01
    } else if(D <= -1.754 | D >= 2.107) {
      sig.level <- 0.05
    } else if(D <= -1.534 | D >= 1.752) {
      sig.level <- 0.1
    }
  } else if(n >= 200) {
    if(D <= -2.462 | D >= 3.492) {
      sig.level <- 0.001
    } else if(D <= -2.132 | D >= 2.768) {
      sig.level <- 0.01
    } else if(D <= -1.760 | D >= 2.100) {
      sig.level <- 0.05
    } else if(D <= -1.539 | D >= 1.748) {
      sig.level <- 0.1
    }
  } else if(n >= 175) {
    if(D <= -2.470 | D >= 3.470) {
      sig.level <- 0.001
    } else if(D <= -2.138 | D >= 2.757) {
      sig.level <- 0.01
    } else if(D <= -1.765 | D >= 2.095) {
      sig.level <- 0.05
    } else if(D <= -1.542 | D >= 1.746) {
      sig.level <- 0.1
    }
  } else if(n >= 150) {
    if(D <= -2.477 | D >= 3.443) {
      sig.level <- 0.001
    } else if(D <= -2.144 | D >= 2.743) {
      sig.level <- 0.01
    } else if(D <= -1.769 | D >= 2.089) {
      sig.level <- 0.05
    } else if(D <= -1.545 | D >= 1.743) {
      sig.level <- 0.1
    }
  } else if(n >= 140) {
    if(D <= -2.481 | D >= 3.430) {
      sig.level <- 0.001
    } else if(D <= -2.147 | D >= 2.736) {
      sig.level <- 0.01
    } else if(D <= -1.771 | D >= 2.086) {
      sig.level <- 0.05
    } else if(D <= -1.547 | D >= 1.741) {
      sig.level <- 0.1
    }
  } else if(n >= 130) {
    if(D <= -2.484 | D >= 3.416) {
      sig.level <- 0.001
    } else if(D <= -2.150 | D >= 2.730) {
      sig.level <- 0.01
    } else if(D <= -1.774 | D >= 2.084) {
      sig.level <- 0.05
    } else if(D <= -1.549 | D >= 1.740) {
      sig.level <- 0.1
    }
  } else if(n >= 120) {
    if(D <= -2.488 | D >= 3.401) {
      sig.level <- 0.001
    } else if(D <= -2.153 | D >= 2.722) {
      sig.level <- 0.01
    } else if(D <= -1.776 | D >= 2.080) {
      sig.level <- 0.05
    } else if(D <= -1.550 | D >= 1.739) {
      sig.level <- 0.1
    }
  } else if(n >= 110) {
    if(D <= -2.492 | D >= 3.385) {
      sig.level <- 0.001
    } else if(D <= -2.157 | D >= 2.713) {
      sig.level <- 0.01
    } else if(D <= -1.779 | D >= 2.077) {
      sig.level <- 0.05
    } else if(D <= -1.552 | D >= 1.737) {
      sig.level <- 0.1
    }
  } else if(n >= 100) {
    if(D <= -2.495 | D >= 3.336) {
      sig.level <- 0.001
    } else if(D <= -2.160 | D >= 2.704) {
      sig.level <- 0.01
    } else if(D <= -1.781 | D >= 2.073) {
      sig.level <- 0.05
    } else if(D <= -1.555 | D >= 1.735) {
      sig.level <- 0.1
    }
  } else if(n >= 95) {
    if(D <= -2.497 | D >= 3.355) {
      sig.level <- 0.001
    } else if(D <= -2.162 | D >= 2.699) {
      sig.level <- 0.01
    } else if(D <= -1.783 | D >= 2.071) {
      sig.level <- 0.05
    } else if(D <= -1.556 | D >= 1.734) {
      sig.level <- 0.1
    }
  } else if(n >= 90) {
    if(D <= -2.499 | D >= 3.345) {
      sig.level <- 0.001
    } else if(D <= -2.164 | D >= 2.693) {
      sig.level <- 0.01
    } else if(D <= -1.784 | D >= 2.069) {
      sig.level <- 0.05
    } else if(D <= -1.557 | D >= 1.733) {
      sig.level <- 0.1
    }
  } else if(n >= 85) {
    if(D <= -2.500 | D >= 3.333) {
      sig.level <- 0.001
    } else if(D <= -2.166 | D >= 2.687) {
      sig.level <- 0.01
    } else if(D <= -1.786 | D >= 2.066) {
      sig.level <- 0.05
    } else if(D <= -1.559 | D >= 1.732) {
      sig.level <- 0.1
    }
  } else if(n >= 80) {
    if(D <= -2.502 | D >= 3.320) {
      sig.level <- 0.001
    } else if(D <= -2.168 | D >= 2.681) {
      sig.level <- 0.01
    } else if(D <= -1.788 | D >= 2.064) {
      sig.level <- 0.05
    } else if(D <= -1.560 | D >= 1.731) {
      sig.level <- 0.1
    }
  } else if(n >= 75) {
    if(D <= -2.504 | D >= 3.306) {
      sig.level <- 0.001
    } else if(D <= -2.170 | D >= 2.673) {
      sig.level <- 0.01
    } else if(D <= -1.790 | D >= 2.061) {
      sig.level <- 0.05
    } else if(D <= -1.561 | D >= 1.730) {
      sig.level <- 0.1
    }
  } else if(n >= 70) {
    if(D <= -2.505 | D >= 3.291) {
      sig.level <- 0.001
    } else if(D <= -2.171 | D >= 2.666) {
      sig.level <- 0.01
    } else if(D <= -1.791 | D >= 2.058) {
      sig.level <- 0.05
    } else if(D <= -1.563 | D >= 1.729) {
      sig.level <- 0.1
    }
  } else if(n >= 65) {
    if(D <= -2.506 | D >= 3.274) {
      sig.level <- 0.001
    } else if(D <= -2.173 | D >= 2.658) {
      sig.level <- 0.01
    } else if(D <= -1.793 | D >= 2.055) {
      sig.level <- 0.05
    } else if(D <= -1.565 | D >= 1.727) {
      sig.level <- 0.1
    }
  } else if(n >= 60) {
    if(D <= -2.506 | D >= 3.256) {
      sig.level <- 0.001
    } else if(D <= -2.175 | D >= 2.649) {
      sig.level <- 0.01
    } else if(D <= -1.795 | D >= 2.052) {
      sig.level <- 0.05
    } else if(D <= -1.566 | D >= 1.726) {
      sig.level <- 0.1
    }
  } else if(n >= 55) {
    if(D <= -2.506 | D >= 3.235) {
      sig.level <- 0.001
    } else if(D <= -2.177 | D >= 2.638) {
      sig.level <- 0.01
    } else if(D <= -1.797 | D >= 2.048) {
      sig.level <- 0.05
    } else if(D <= -1.568 | D >= 1.724) {
      sig.level <- 0.1
    }
  } else if(n >= 50) {
    if(D <= -2.505 | D >= 3.212) {
      sig.level <- 0.001
    } else if(D <= -2.178 | D >= 2.627) {
      sig.level <- 0.01
    } else if(D <= -1.800 | D >= 2.044) {
      sig.level <- 0.05
    } else if(D <= -1.570 | D >= 1.723) {
      sig.level <- 0.1
    }
  } else if(n >= 49) {
    if(D <= -2.505 | D >= 3.207) {
      sig.level <- 0.001
    } else if(D <= -2.178 | D >= 2.624) {
      sig.level <- 0.01
    } else if(D <= -1.800 | D >= 2.042) {
      sig.level <- 0.05
    } else if(D <= -1.571 | D >= 1.722) {
      sig.level <- 0.1
    }
  } else if(n >= 48) {
    if(D <= -2.505 | D >= 3.202) {
      sig.level <- 0.001
    } else if(D <= -2.178 | D >= 2.622) {
      sig.level <- 0.01
    } else if(D <= -1.800 | D >= 2.042) {
      sig.level <- 0.05
    } else if(D <= -1.571 | D >= 1.722) {
      sig.level <- 0.1
    }
  } else if(n >= 47) {
    if(D <= -2.504 | D >= 3.196) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.619) {
      sig.level <- 0.01
    } else if(D <= -1.801 | D >= 2.041) {
      sig.level <- 0.05
    } else if(D <= -1.572 | D >= 1.722) {
      sig.level <- 0.1
    }
  } else if(n >= 46) {
    if(D <= -2.504 | D >= 3.191) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.617) {
      sig.level <- 0.01
    } else if(D <= -1.801 | D >= 2.040) {
      sig.level <- 0.05
    } else if(D <= -1.572 | D >= 1.721) {
      sig.level <- 0.1
    }
  } else if(n >= 45) {
    if(D <= -2.503 | D >= 3.185) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.613) {
      sig.level <- 0.01
    } else if(D <= -1.802 | D >= 2.039) {
      sig.level <- 0.05
    } else if(D <= -1.573 | D >= 1.721) {
      sig.level <- 0.1
    }
  } else if(n >= 44) {
    if(D <= -2.502 | D >= 3.180) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.611) {
      sig.level <- 0.01
    } else if(D <= -1.802 | D >= 2.038) {
      sig.level <- 0.05
    } else if(D <= -1.573 | D >= 1.721) {
      sig.level <- 0.1
    }
  } else if(n >= 43) {
    if(D <= -2.502 | D >= 3.173) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.608) {
      sig.level <- 0.01
    } else if(D <= -1.803 | D >= 2.037) {
      sig.level <- 0.05
    } else if(D <= -1.574 | D >= 1.720) {
      sig.level <- 0.1
    }
  } else if(n >= 42) {
    if(D <= -2.501 | D >= 3.168) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.605) {
      sig.level <- 0.01
    } else if(D <= -1.803 | D >= 2.036) {
      sig.level <- 0.05
    } else if(D <= -1.574 | D >= 1.720) {
      sig.level <- 0.1
    }
  } else if(n >= 41) {
    if(D <= -2.500 | D >= 3.160) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.601) {
      sig.level <- 0.01
    } else if(D <= -1.803 | D >= 2.034) {
      sig.level <- 0.05
    } else if(D <= -1.575 | D >= 1.719) {
      sig.level <- 0.1
    }
  } else if(n >= 40) {
    if(D <= -2.499 | D >= 3.155) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.599) {
      sig.level <- 0.01
    } else if(D <= -1.804 | D >= 2.033) {
      sig.level <- 0.05
    } else if(D <= -1.575 | D >= 1.719) {
      sig.level <- 0.1
    }
  } else if(n >= 39) {
    if(D <= -2.498 | D >= 3.147) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.595) {
      sig.level <- 0.01
    } else if(D <= -1.804 | D >= 2.032) {
      sig.level <- 0.05
    } else if(D <= -1.576 | D >= 1.718) {
      sig.level <- 0.1
    }
  } else if(n >= 38) {
    if(D <= -2.496 | D >= 3.141) {
      sig.level <- 0.001
    } else if(D <= -2.178 | D >= 2.592) {
      sig.level <- 0.01
    } else if(D <= -1.804 | D >= 2.031) {
      sig.level <- 0.05
    } else if(D <= -1.576 | D >= 1.718) {
      sig.level <- 0.1
    }
  } else if(n >= 37) {
    if(D <= -2.495 | D >= 3.132) {
      sig.level <- 0.001
    } else if(D <= -2.179 | D >= 2.588) {
      sig.level <- 0.01
    } else if(D <= -1.805 | D >= 2.030) {
      sig.level <- 0.05
    } else if(D <= -1.577 | D >= 1.717) {
      sig.level <- 0.1
    }
  } else if(n >= 36) {
    if(D <= -2.493 | D >= 3.126) {
      sig.level <- 0.001
    } else if(D <= -2.178 | D >= 2.585) {
      sig.level <- 0.01
    } else if(D <= -1.805 | D >= 2.029) {
      sig.level <- 0.05
    } else if(D <= -1.577 | D >= 1.717) {
      sig.level <- 0.1
    }
  } else if(n >= 35) {
    if(D <= -2.492 | D >= 3.116) {
      sig.level <- 0.001
    } else if(D <= -2.178 | D >= 2.580) {
      sig.level <- 0.01
    } else if(D <= -1.806 | D >= 2.027) {
      sig.level <- 0.05
    } else if(D <= -1.578 | D >= 1.717) {
      sig.level <- 0.1
    }
  } else if(n >= 34) {
    if(D <= -2.489 | D >= 3.110) {
      sig.level <- 0.001
    } else if(D <= -2.177 | D >= 2.577) {
      sig.level <- 0.01
    } else if(D <= -1.806 | D >= 2.026) {
      sig.level <- 0.05
    } else if(D <= -1.578 | D >= 1.716) {
      sig.level <- 0.1
    }
  } else if(n >= 33) {
    if(D <= -2.487 | D >= 3.099) {
      sig.level <- 0.001
    } else if(D <= -2.177 | D >= 2.572) {
      sig.level <- 0.01
    } else if(D <= -1.806 | D >= 2.024) {
      sig.level <- 0.05
    } else if(D <= -1.579 | D >= 1.716) {
      sig.level <- 0.1
    }
  } else if(n >= 32) {
    if(D <= -2.484 | D >= 3.092) {
      sig.level <- 0.001
    } else if(D <= -2.175 | D >= 2.569) {
      sig.level <- 0.01
    } else if(D <= -1.806 | D >= 2.023) {
      sig.level <- 0.05
    } else if(D <= -1.579 | D >= 1.715) {
      sig.level <- 0.1
    }
  } else if(n >= 31) {
    if(D <= -2.482 | D >= 3.080) {
      sig.level <- 0.001
    } else if(D <= -2.175 | D >= 2.563) {
      sig.level <- 0.01
    } else if(D <= -1.807 | D >= 2.021) {
      sig.level <- 0.05
    } else if(D <= -1.580 | D >= 1.714) {
      sig.level <- 0.1
    }
  } else if(n >= 30) {
    if(D <= -2.478 | D >= 3.073) {
      sig.level <- 0.001
    } else if(D <= -2.173 | D >= 2.559) {
      sig.level <- 0.01
    } else if(D <= -1.807 | D >= 2.020) {
      sig.level <- 0.05
    } else if(D <= -1.580 | D >= 1.714) {
      sig.level <- 0.1
    }
  } else if(n >= 29) {
    if(D <= -2.475 | D >= 3.060) {
      sig.level <- 0.001
    } else if(D <= -2.173 | D >= 2.553) {
      sig.level <- 0.01
    } else if(D <= -1.807 | D >= 2.018) {
      sig.level <- 0.05
    } else if(D <= -1.581 | D >= 1.714) {
      sig.level <- 0.1
    }
  } else if(n >= 28) {
    if(D <= -2.471 | D >= 3.052) {
      sig.level <- 0.001
    } else if(D <= -2.171 | D >= 2.549) {
      sig.level <- 0.01
    } else if(D <= -1.807 | D >= 2.017) {
      sig.level <- 0.05
    } else if(D <= -1.581 | D >= 1.713) {
      sig.level <- 0.1
    }
  } else if(n >= 27) {
    if(D <= -2.467 | D >= 3.037) {
      sig.level <- 0.001
    } else if(D <= -2.170 | D >= 2.542) {
      sig.level <- 0.01
    } else if(D <= -1.807 | D >= 2.014) {
      sig.level <- 0.05
    } else if(D <= -1.582 | D >= 1.712) {
      sig.level <- 0.1
    }
  } else if(n >= 26) {
    if(D <= -2.461 | D >= 3.029) {
      sig.level <- 0.001
    } else if(D <= -2.167 | D >= 2.538) {
      sig.level <- 0.01
    } else if(D <= -1.807 | D >= 2.013) {
      sig.level <- 0.05
    } else if(D <= -1.582 | D >= 1.712) {
      sig.level <- 0.1
    }
  } else if(n >= 25) {
    if(D <= -2.457 | D >= 3.011) {
      sig.level <- 0.001
    } else if(D <= -2.165 | D >= 2.530) {
      sig.level <- 0.01
    } else if(D <= -1.807 | D >= 2.010) {
      sig.level <- 0.05
    } else if(D <= -1.583 | D >= 1.712) {
      sig.level <- 0.1
    }
  } else if(n >= 24) {
    if(D <= -2.449 | D >= 3.002) {
      sig.level <- 0.001
    } else if(D <= -2.162 | D >= 2.526) {
      sig.level <- 0.01
    } else if(D <= -1.806 | D >= 2.009) {
      sig.level <- 0.05
    } else if(D <= -1.583 | D >= 1.712) {
      sig.level <- 0.1
    }
  } else if(n >= 23) {
    if(D <= -2.443 | D >= 2.983) {
      sig.level <- 0.001
    } else if(D <= -2.160 | D >= 2.516) {
      sig.level <- 0.01
    } else if(D <= -1.806 | D >= 2.006) {
      sig.level <- 0.05
    } else if(D <= -1.584 | D >= 1.710) {
      sig.level <- 0.1
    }
  } else if(n >= 22) {
    if(D <= -2.434 | D >= 2.973) {
      sig.level <- 0.001
    } else if(D <= -2.153 | D >= 2.512) {
      sig.level <- 0.01
    } else if(D <= -1.804 | D >= 2.005) {
      sig.level <- 0.05
    } else if(D <= -1.584 | D >= 1.711) {
      sig.level <- 0.1
    }
  } else if(n >= 21) {
    if(D <= -2.426 | D >= 2.950) {
      sig.level <- 0.001
    } else if(D <= -2.152 | D >= 2.501) {
      sig.level <- 0.01
    } else if(D <= -1.805 | D >= 2.001) {
      sig.level <- 0.05
    } else if(D <= -1.585 | D >= 1.709) {
      sig.level <- 0.1
    }
  } else if(n >= 20) {
    if(D <= -2.414 | D >= 2.939) {
      sig.level <- 0.001
    } else if(D <= -2.146 | D >= 2.496) {
      sig.level <- 0.01
    } else if(D <= -1.803 | D >= 2.001) {
      sig.level <- 0.05
    } else if(D <= -1.584 | D >= 1.710) {
      sig.level <- 0.1
    }
  } else if(n >= 19) {
    if(D <= -2.403 | D >= 2.911) {
      sig.level <- 0.001
    } else if(D <= -2.141 | D >= 2.483) {
      sig.level <- 0.01
    } else if(D <= -1.802 | D >= 1.996) {
      sig.level <- 0.05
    } else if(D <= -1.585 | D >= 1.708) {
      sig.level <- 0.1
    }
  } else if(n >= 18) {
    if(D <= -2.387 | D >= 2.900) {
      sig.level <- 0.001
    } else if(D <= -2.132 | D >= 2.478) {
      sig.level <- 0.01
    } else if(D <= -1.799 | D >= 1.996) {
      sig.level <- 0.05
    } else if(D <= -1.584 | D >= 1.709) {
      sig.level <- 0.1
    }
  } else if(n >= 17) {
    if(D <= -2.372 | D >= 2.866) {
      sig.level <- 0.001
    } else if(D <= -2.126 | D >= 2.461) {
      sig.level <- 0.01
    } else if(D <= -1.798 | D >= 1.990) {
      sig.level <- 0.05
    } else if(D <= -1.585 | D >= 1.708) {
      sig.level <- 0.1
    }
  } else if(n >= 16) {
    if(D <= -2.350 | D >= 2.854) {
      sig.level <- 0.001
    } else if(D <= -2.113 | D >= 2.457) {
      sig.level <- 0.01
    } else if(D <= -1.793 | D >= 1.990) {
      sig.level <- 0.05
    } else if(D <= -1.583 | D >= 1.709) {
      sig.level <- 0.1
    }
  } else if(n >= 15) {
    if(D <= -2.329 | D >= 2.811) {
      sig.level <- 0.001
    } else if(D <= -2.103 | D >= 2.436) {
      sig.level <- 0.01
    } else if(D <= -1.791 | D >= 1.984) {
      sig.level <- 0.05
    } else if(D <= -1.584 | D >= 1.708) {
      sig.level <- 0.1
    }
  } else if(n >= 14) {
    if(D <= -2.299 | D >= 2.798) {
      sig.level <- 0.001
    } else if(D <= -2.085 | D >= 2.432) {
      sig.level <- 0.01
    } else if(D <= -1.783 | D >= 1.985) {
      sig.level <- 0.05
    } else if(D <= -1.580 | D >= 1.710) {
      sig.level <- 0.1
    }
  } else if(n >= 13) {
    if(D <= -2.267 | D >= 2.741) {
      sig.level <- 0.001
    } else if(D <= -2.069 | D >= 2.403) {
      sig.level <- 0.01
    } else if(D <= -1.779 | D >= 1.976) {
      sig.level <- 0.05
    } else if(D <= -1.580 | D >= 1.708) {
      sig.level <- 0.1
    }
  } else if(n >= 12) {
    if(D <= -2.223 | D >= 2.729) {
      sig.level <- 0.001
    } else if(D <= -2.041 | D >= 2.401) {
      sig.level <- 0.01
    } else if(D <= -1.765 | D >= 1.979) {
      sig.level <- 0.05
    } else if(D <= -1.573 | D >= 1.713) {
      sig.level <- 0.1
    }
  } else if(n >= 11) {
    if(D <= -2.174 | D >= 2.649) {
      sig.level <- 0.001
    } else if(D <= -2.014 | D >= 2.359) {
      sig.level <- 0.01
    } else if(D <= -1.757 | D >= 1.966) {
      sig.level <- 0.05
    } else if(D <= -1.572 | D >= 1.710) {
      sig.level <- 0.1
    }
  } else if(n >= 10) {
    if(D <= -2.105 | D >= 2.640) {
      sig.level <- 0.001
    } else if(D <= -1.967 | D >= 2.362) {
      sig.level <- 0.01
    } else if(D <= -1.733 | D >= 1.975) {
      sig.level <- 0.05
    } else if(D <= -1.559 | D >= 1.719) {
      sig.level <- 0.1
    }
  } else if(n >= 9) {
    if(D <= -2.023 | D >= 2.519) {
      sig.level <- 0.001
    } else if(D <= -1.916 | D >= 2.296) {
      sig.level <- 0.01
    } else if(D <= -1.713 | D >= 1.954) {
      sig.level <- 0.05
    } else if(D <= -1.553 | D >= 1.715) {
      sig.level <- 0.1
    }
  } else if(n >= 8) {
    if(D <= -1.909 | D >= 2.524) {
      sig.level <- 0.001
    } else if(D <= -1.830 | D >= 2.313) {
      sig.level <- 0.01
    } else if(D <= -1.663 | D >= 1.975) {
      sig.level <- 0.05
    } else if(D <= -1.522 | D >= 1.736) {
      sig.level <- 0.1
    }
  } else if(n >= 7) {
    if(D <= -1.761 | D >= 2.524) {
      sig.level <- 0.001
    } else if(D <= -1.721 | D >= 2.185) {
      sig.level <- 0.01
    } else if(D <= -1.608 | D >= 1.932) {
      sig.level <- 0.05
    } else if(D <= -1.498 | D >= 1.728) {
      sig.level <- 0.1
    }
  } else if(n >= 6) {
    if(D <= -1.556 | D >= 2.373) {
      sig.level <- 0.001
    } else if(D <= -1.540 | D >= 2.255) {
      sig.level <- 0.01
    } else if(D <= -1.478 | D >= 1.999) {
      sig.level <- 0.05
    } else if(D <= -1.405 | D >= 1.786) {
      sig.level <- 0.1
    }
  } else if(n >= 5) {
    if(D <= -1.276 | D >= 1.913) {
      sig.level <- 0.001
    } else if(D <= -1.275 | D >= 1.901) {
      sig.level <- 0.01
    } else if(D <= -1.269 | D >= 1.834) {
      sig.level <- 0.05
    } else if(D <= -1.255 | D >= 1.737) {
      sig.level <- 0.1
    }
  } else if(n == 4) {
    if(D <= -0.876 | D >= 2.336) {
      sig.level <- 0.001
    } else if(D <= -0.876 | D >= 2.324) {
      sig.level <- 0.01
    } else if(D <= -0.876 | D >= 2.232) {
      sig.level <- 0.05
    } else if(D <= -0.876 | D >= 2.081) {
      sig.level <- 0.1
    }
  } 
  
  #return(c(D, sig.level))
  cat("D\tsig.level\n", sep="")
  cat(D, "\t", sig.level, "\n", sep="")
  #print(D, sig.level)
  
}


