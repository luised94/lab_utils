# X509_V__FLAG_CB_ISSUER_CHECK

> Tried to download tools to handle pdfs (pdf2doi,pdf2txt,pdfrenamer,pdf2bib) and ran into the following error:

AttributeError: module 'lib' has no attribute 'X509_V_FLAG_CB_ISSUER_CHECK'

> Conflict between openssl and cryptography program

> See https://github.com/pyca/cryptography/issues/7126 and https://github.com/Azure/azure-cli/issues/16858

> Ran the following:

'$python3 -m pip install -U pip '

> Was able to install pdf2doi, pdf2bib, pdf-renamer afterwards

'$pip install pdf-renamer==1.0rc6 '

> On my second computer, had to install/update pip manually the following:

'$curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py'
'$python3 get-pip.py'