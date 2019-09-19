.. _citation:

-------------
Citing abagen
-------------

We're thrilled you've found ``abagen`` useful in your work! To make citing
it as easy as possible we have created a `DOI`_ with `Zenodo`_. Since much of
the current codebase was inspired by the workflow laid out in Arnatkevičiūte et
al., 2019, and it all uses data released with Hawrylycz et al., 2012, we ask
that you please cite those papers as well:

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
         getLatestIDFromconceptID("", cb);
      } else {
         getZenodoIDFromTag("", version, cb);
      }
   }
   </script>

   <p>
     1. <span id="abagen_citation">abagen</span> Available from: <a id="abagen_doi_url" href="https://doi.org/xx.xxxx/zenodo.xxxxxx">xx.xxxx/zenodo.xxxxxx</a>.
     <img src onerror='fillCitation()' alt="" />
   </p>
   <p>
     2. Arnatkevic̆iūtė, A., Fulcher, B. D., & Fornito, A. (2019). A practical guide to linking brain-wide gene expression and neuroimaging data. NeuroImage, 189, 353-367. doi:10.1016/j.neuroimage.2019.01.011
   </p>
   <p>
     3. Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng, L., Miller, J. A., … Jones, A. R. (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489(7416), 391–399. doi:10.1038/nature11405
   </p>

.. _DOI: https://en.wikipedia.org/wiki/Digital_object_identifier
.. _Zenodo: https://zenodo.org
