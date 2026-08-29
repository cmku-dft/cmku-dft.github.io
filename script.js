const lightbox = document.querySelector('#activityLightbox');

if (lightbox) {
  const enlargedPhoto = lightbox.querySelector('img');
  const closeButton = lightbox.querySelector('.lightboxClose');

  document.querySelectorAll('.galleryButton').forEach((button) => {
    button.addEventListener('click', () => {
      const thumbnail = button.querySelector('img');
      enlargedPhoto.src = button.dataset.full;
      enlargedPhoto.alt = thumbnail.alt;
      enlargedPhoto.classList.toggle('portrait', thumbnail.naturalHeight > thumbnail.naturalWidth);
      lightbox.showModal();
    });
  });

  closeButton.addEventListener('click', () => lightbox.close());
  lightbox.addEventListener('click', (event) => {
    if (event.target === lightbox) lightbox.close();
  });
}
