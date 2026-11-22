## Overview

We propose an algorithm that automatically identifies the tree topology of an electrical network from telemetry data, using current or power measurements collected by sensors installed on a subset of its feeders.

The method first exploits the ordering of current magnitudes to determine the set of candidate parents for each monitored feeder. It then selects the most likely parent by maximizing the correlation between current variations.

The parent-selection problem is formulated as an integer program. Its linear programming (LP) relaxation provides an upper bound on the optimal objective value and, for the base problem, is tight. We also introduce an extension incorporating capacity constraints.
