"""Verified HTTPS, including the Windows system trust store used by this workstation."""
import ssl
import urllib.request


def open_url(request, timeout=120):
    context = ssl.create_default_context()
    if hasattr(ssl, 'enum_certificates'):
        for store in ('ROOT', 'CA'):
            for cert, encoding, trust in ssl.enum_certificates(store):
                if encoding == 'x509_asn' and (trust is True or ssl.Purpose.SERVER_AUTH.oid in trust):
                    context.load_verify_locations(cadata=ssl.DER_cert_to_PEM_cert(cert))
    return urllib.request.urlopen(request, timeout=timeout, context=context)
