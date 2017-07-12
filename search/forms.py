from django import forms


# A form where users can contact the website Admin
class ContactForm(forms.Form):
    name = forms.CharField(required=True, max_length=100, help_text='100 chars max')
    email = forms.EmailField(required=True)
    comment = forms.CharField(required=True, widget=forms.Textarea)
