import requests


def get_authentication_token(nomad_url, username, password):
    '''Get the token for accessing your NOMAD unpublished uploads remotely'''
    try:
        response = requests.get(
            nomad_url + 'auth/token', params=dict(username=username, password=password), timeout=10)
        token = response.json().get('access_token')
        if token:
            return token

        raise ValueError('Noda response is missing token: \n' + response.json)
    except Exception:
        print('something went wrong trying to get authentication token')
        return


def create_dataset(nomad_url, token, dataset_name):
    '''Create a dataset to group a series of NOMAD entries'''
    try:
        response = requests.post(
            nomad_url + 'datasets/',
            headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
            json={"dataset_name": dataset_name},
            timeout=10
        )
        dataset_id = response.json().get('dataset_id')
        if dataset_id:
            return dataset_id

        raise ValueError('response is missing dataset_id:\n ' + response.json())
        return
    except Exception:
        print('something went wrong trying to create a dataset')
        return


def upload_to_nomad(nomad_url, token, upload_file):
    '''Upload a single file for NOMAD upload, e.g., zip format'''
    with open(upload_file, 'rb') as f:
        try:
            response = requests.post(
                nomad_url + 'uploads',
                headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
                data=f, timeout=30)
            upload_id = response.json().get('upload_id')
            if upload_id:
                return upload_id
            raise ValueError('response is missing upload_id:\n ' + response.json())
            return
        except Exception:
            print('something went wrong uploading to NOMAD')
            return


def check_upload_status(nomad_url, token, upload_id):
    '''
    # upload success => returns 'Process publish_upload completed successfully'
    # publish success => 'Process publish_upload completed successfully'
    '''
    try:
        response = requests.get(
            nomad_url + 'uploads/' + upload_id,
            headers={'Authorization': f'Bearer {token}'}, timeout=30)
        status_message = response.json().get('data').get('last_status_message')
        if status_message:
            return status_message

        print('response is missing status_message: ')
        print(response.json())
        return
    except Exception:
        print('something went wrong trying to check the status of upload' + upload_id)
        # upload gets deleted from the upload staging area once published...or in this case something went wrong
        return


def edit_upload_metadata(nomad_url, token, upload_id, metadata):
    '''
    Example of new metadata:
    upload_name = 'Test_Upload_Name'
    metadata = {
        "metadata": {
        "upload_name": upload_name,
        "references": ["https://doi.org/xx.xxxx/xxxxxx"],
        "datasets": dataset_id,
        "embargo_length": 0,
        "coauthors": ["coauthor@affiliation.de"],
        "comment": 'This is a test upload...'
        },
    }
    '''

    try:
        response = requests.post(
            nomad_url + 'uploads/' + upload_id + '/edit',
            headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
            json=metadata, timeout=30)
        return response
    except Exception:
        print('something went wrong trying to add metadata to upload' + upload_id)
        return


def publish_upload(nomad_url, token, upload_id):
    '''Publish an upload'''
    try:
        response = requests.post(
            nomad_url + 'uploads/' + upload_id + '/action/publish',
            headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
            timeout=30)
        return response
    except Exception:
        print('something went wrong trying to publish upload: ' + upload_id)
        return
