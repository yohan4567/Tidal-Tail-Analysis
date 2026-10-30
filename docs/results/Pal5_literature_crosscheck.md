# Pal 5: literature-rate cross-check

The current conditional grid does not establish a physical rate of 31 solar
masses/Myr. That is the best sampled candidate and lies on the rate upper edge.
The grid includes a rate of 5; it has profile q=2.921929361823686, with its best
sampled time 1,100 Myr and escape-mass tilt +1. The global grid best is rate 31,
time 900 Myr and tilt -1. This demonstrates sensitivity to the assumed escaping
population. No calibrated confidence set or goodness-of-fit conclusion follows
from these numbers alone, and the known literature estimate is not imposed as a
fitting target or used to discard unfavourable models.

Kuepper et al. (2015), Section IV.2.4, summarise Odenkirchen et al. (2003): using
an assumed present-day mass of 6,000 solar masses, estimates from trailing and
leading star counts were 4.9 and 3.8 solar masses/Myr, respectively. Their own
weakly constrained streakline inference is 7.9 (+5.3/-6.4) solar masses/Myr.
They also quote a fractional loss estimate of (0.7 +/-0.2)/Gyr: at their adopted
16,000 solar masses this scales to 11.2 +/-3.2 solar masses/Myr. These quantities
are model-dependent estimates rather than a universal exact 5 solar masses/Myr.

Sources:
- Kuepper et al. 2015, Section IV.2.4: https://arxiv.org/html/1502.02658
- Odenkirchen et al. 2003: https://arxiv.org/abs/astro-ph/0307446

Our adopted present mass is 13,400 solar masses; our rate counts the mass
represented by living-star population nodes, assumes f=1 for the prescribed
history and fits Gaia-selected counts in a masked field. Literature comparisons
need matched mass basis, completeness, stellar mass functions, background,
field/time coverage and dynamical assumptions. The population tilt changes the
conversion from selected counts to escaping mass; the halo/orbit and time history
are further systematic uncertainties. Before reporting an absolute rate, run
population/CMD and selection-transfer checks, spatial/background diagnostics,
numerical refinements and calibrated held-out recovery. Do not simply enlarge the
upper boundary until an interior best point appears, or force agreement with 5.
