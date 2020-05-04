docker run \
       -it \
       matthewvc1/primer_eval:0.1.2 \
       -v databases:/usr/local/databases \
       -v output:/home/runDir/output \
       --name primerEval \
       bash
