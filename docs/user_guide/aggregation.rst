.. _usage_aggregation:

Sample aggregation options
==========================

The primary goal of :func:`abagen.get_expression_data` is to allow users to
aggregate the ~3,700 disparate tissue samples from the Allen Human Brain Atlas
into regions of interest defined by an atlas or parcellation file. However,
there exist several options for exactly *how* to aggregate samples within each
region of the specified atlas.

These options are controlled via two parameters to
:func:`abagen.get_expression_data`: ``region_agg`` and ``agg_metric``. We
discuss both parameters and the different options available to each below.

.. _usage_region_agg:

The `region_agg` parameter
--------------------------

This parameter determines how samples are aggregated together to generate the
expression values for a region. It can take two values: `'donors'` or
`'samples'`.

If set to `'donors'`, expression values for all samples assigned to a region
are aggregated independently for each donor and *then* aggregated across
donors. If set to `'samples'`, expression values for all samples for all
donors assigned to a region are aggregated simultaneously.

.. _usage_agg_metric:

The `agg_metric` parameter
--------------------------

This parameter determines the actual metric used for aggregating samples into
regional expression values. It can be set to any callable function (as long as
that function accepts the keyword ``axis`` argument), but generally either
`'mean'` (the default) or `'median'` will suffice.
