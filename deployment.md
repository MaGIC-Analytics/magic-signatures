# Deployment notes — magic-signature
Based on rand3k repo

The normal path is **CI/CD**: Cloud Run is connected to the `MaGIC-Analytics/magic-signature`
GitHub repo and builds + deploys on push (the trigger lives in GCP). The CLI/manual steps below
are a fallback. Note the scoring engines (GSVA/AUCell/singscore) are CPU/RAM-heavy — allocate
generous memory (≥ 2–4 GB) and a longer build timeout for the Bioconductor install.

## CLI deployment
```
PROJECTID=$(gcloud config get-value project)
docker build . -t gcr.io/$PROJECTID/magic-signature
docker push gcr.io/$PROJECTID/magic-signature
gcloud run deploy --image gcr.io/$PROJECTID/magic-signature --platform managed --max-instances 1
```
Manually adjust CPUs and RAM applied to the container as it may be custom.


## Manual deployment
Also very easy:
- Log into your gcloud account. 
- Access cloudrun
- Hit deploy new.
- Choose the repository and set it to the dockerfile
- Define your build conditions
- Update the timeout. Default is 10 min. 
- Once built, add mapping to cloud run to add that tool
    - Copy the cname, go to cloud dns under network services. Add the record set for the new tool