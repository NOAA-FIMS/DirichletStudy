# Hake Functional Analysis Results

This analysis compares the linear and saturated Dirichlet variants on the generated hake composition output. The input data aggregate to counts of 75,778, 75,807, and 74,283 across the three composition categories, giving observed proportions of 0.335497, 0.335625, and 0.328878. The functional analysis evaluated 231 open-simplex points generated from `h = 0.05`.

## Main Result

Both variants identify the same best-supported composition on the evaluated grid:

| variant | best p1 | best p2 | best p3 | max loglik | L1 distance | L2 distance |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| linear | 0.347826 | 0.347826 | 0.304348 | -170.279 | 0.049060 | 0.030043 |
| saturated | 0.347826 | 0.347826 | 0.304348 | -25.435 | 0.049060 | 0.030043 |

The selected grid point is close to the observed composition, especially given the simplex spacing. This is useful for a fisheries modeling argument: the analysis recovers the empirical age or length composition pattern without requiring the optimum to sit exactly on the observed proportions.

## Linear Variant

The linear variant produces a sharply peaked likelihood surface. Its maximum log-likelihood is -170.279, while the minimum is -89,992.028. The worst grid point is near the boundary at `(0.043478, 0.043478, 0.913043)`, with L1 distance 1.168331 and L2 distance 0.715454 from the observed composition.

This strong contrast means the linear form is highly discriminating. It heavily penalizes compositions that place most mass in one category when the hake output is nearly balanced across categories. In a stock-assessment setting, that behavior is attractive when the goal is strong information content from composition data, but it can also imply sensitivity to aggregation, effective sample size, and model misspecification.

The linear effective sample size diagnostic is 112,934.5. That scale is consistent with a very concentrated likelihood surface and very large gradients. The mean gradient norm is 403,005.9, and the mean absolute derivative step is 40,005.0 in the summary file. Directional derivative diagnostics are also large, especially for the first two coordinates, indicating that small moves across the simplex can create large likelihood changes.

## Saturated Variant

The saturated variant gives the same best grid point but a much flatter likelihood surface. Its maximum log-likelihood is -25.435 and its minimum is -30.225. The worst grid point is again `(0.043478, 0.043478, 0.913043)`, but its relative likelihood remains 0.008314 rather than collapsing to zero.

This makes the saturated form a compelling robustness contrast for fisheries applications. It still favors the composition nearest the hake output, but it avoids assigning overwhelming penalty to far-away compositions. That can be useful when composition samples are overdispersed, have unmodeled haul effects, or include process variation not captured by a simpler multinomial-like likelihood.

The saturated effective sample size diagnostic is approximately 2.0, which matches the flatter likelihood behavior. The mean gradient norm is 39.1, and the mean absolute derivative step is 5.39. The derivative total variation is much smaller than in the linear variant, suggesting a smoother and less brittle response over the same simplex grid.

## Boundary Behavior

The boundary gradient ratio is 1.56 for the linear variant and 2.23 for the saturated variant, using a boundary threshold of 0.05. Both variants show stronger gradients near the simplex boundary than in the interior. That is expected for compositional likelihoods, but it is still worth reporting because boundary behavior matters when observed or predicted fishery compositions contain rare age or length classes.

The saturated model has the larger boundary gradient ratio even though its absolute gradients are much smaller. Interpreted carefully, this means saturation compresses the overall likelihood scale while still preserving stronger relative sensitivity near low-proportion categories.

## Paper-Oriented Interpretation

These results support a clear comparison:

- The linear variant behaves as an information-rich composition likelihood, sharply identifying the observed hake composition and strongly penalizing implausible simplex regions.
- The saturated variant preserves the same modal composition but substantially flattens the likelihood surface, making it a useful robustness or overdispersion-oriented alternative.
- The shared optimum suggests the saturated variant is not moving the fitted composition away from the empirical signal; instead, it changes the strength of evidence assigned to departures from that signal.
- Boundary-gradient diagnostics are relevant for fisheries modeling because rare composition bins often drive instability in likelihood-based estimation.

For a manuscript, the most compelling framing is that the saturated variant offers a controlled way to reduce excessive composition-data influence while retaining the same qualitative information about the best-supported composition. The linear variant can serve as the sharper baseline, while the saturated variant demonstrates robustness to the kinds of overdispersion and boundary sensitivity that commonly appear in fisheries age- and length-composition data.

## Caveat

The optimum is grid-based, not a continuous optimizer result. Because `h = 0.05`, the nearest evaluated point to the observed composition is approximate. A finer simplex or a continuous optimization pass would be useful if the goal is precise parameter estimation rather than functional comparison of likelihood behavior.
