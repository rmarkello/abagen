.. _citation:

-------------
Citing abagen
-------------

.. note::

    We strongly encourage you to use the :ref:`automatically-generated methods
    reports <usage_reporting>` built into ``abagen``, which will provide a list
    of citations based on your selected processing options!

We're thrilled you've found ``abagen`` useful in your work! Please cite the
following manuscripts when referencing your use of the toolbox:

1. Markello, RD, Arnatkevic̆iūtė, A, Poline, J-B, Fulcher, BD, Fornito, A, &
   Misic, B. (2021). Standardizing workflows in imaging transcriptomics with
   the abagen toolbox. Biorxiv. doi:`10.1101/2021.07.08.451635 <https://doi.org/
   10.1101/2021.07.08.451635>`__
2. Arnatkevic̆iūtė, A, Fulcher, BD, & Fornito, A. (2019). A practical guide to
   linking brain-wide gene expression and neuroimaging data. NeuroImage, 189,
   353-367. doi:`10.1016/j.neuroimage.2019.01.011 <https://doi.org/10.1016/
   j.neuroimage.2019.01.011>`__
3. Hawrylycz, MJ, Lein, ES, Guillozet-Bongaarts, AL, Shen, EH, Ng, L, Miller,
   JA, …, & Jones, AR. (2012). An anatomically comprehensive atlas of the adult
   human brain transcriptome. Nature, 489(7416), 391–399.
   doi:`10.1038/nature11405 <https://doi.org/10.1038/nature11405>`__

Additionally, to cite the specific version of the toolbox used in your analyses
you can use the following Zenodo reference:

.. raw:: html

    <script language="javascript">
    var version = 'latest';
    function fillCitation(){
       $('#abagen_version').text(version);

       function cb(err, zenodoID) {
          getCitation(zenodoID, 'vancouver-brackets-no-et-al', function(err, citation) {
             $('#abagen_citation').text(citation);
          });
          getDOI(zenodoID, function(err, DOI) {
             $('#abagen_doi_url').text('https://doi.org/' + DOI);
             $('#abagen_doi_url').attr('href', 'https://doi.org/' + DOI);
          });
       }

       if(version == 'latest') {
          getLatestIDFromconceptID("3451463", cb);
       } else {
          getZenodoIDFromTag("3451463", version, cb);
       }
    }
    </script>

    <p style="margin-left: 30px">
      <span id="abagen_citation">abagen</span> Available from: <a id="abagen_doi_url" href="https://doi.org/10.5281/zenodo.3451463">10.5281/zenodo.3451463</a>.
      <img src onerror='fillCitation()' alt="" />
    </p>

Note that this will always point to the most recent ``abagen`` release; for
older releases please refer to the `Zenodo listing <https://zenodo.org/search?
page=1&size=20&q=conceptrecid:%223451463%22&sort=-version&all_versions=True>`__.

For more information about why citing software is important please refer to
`this article <https://www.software.ac.uk/how-cite-software>`_ from the
Software Sustainability Institute.

.. _DOI: https://en.wikipedia.org/wiki/Digital_object_identifier
.. _Zenodo: https://zenodo.org
